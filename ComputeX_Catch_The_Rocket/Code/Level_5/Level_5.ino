#include <Arduino.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define TRACKER_GLITCHES 1

// ============================================================
// CONFIGURATION
// ============================================================

static const uint32_t TOTAL_READINGS = 100000;
static const uint32_t BLOCK_LEN      = 1000;
static const uint32_t NUM_BLOCKS     = TOTAL_READINGS / BLOCK_LEN;

static const uint32_t SURE_PARACHUTE_BLOCK = 45;
static const int32_t GLITCH_JUMP_CM = 2500;

// ============================================================
// MATHEMATICAL CONSTANTS
// ============================================================

static const double BINOM[5][5] = {
  {1, 0, 0, 0, 0},
  {1, 1, 0, 0, 0},
  {1, 2, 1, 0, 0},
  {1, 3, 3, 1, 0},
  {1, 4, 6, 4, 1}
};

static const double K1000[5] = {
  1.0,
  1e3,
  1e6,
  1e9,
  1e12
};

// ============================================================
// DATA STRUCTURES
// ============================================================

typedef struct {
  int64_t s0;
  int64_t s1;
  int64_t s2;
} BlockSums;

typedef struct {
  double c0;
  double c1;
  double c2;
  double xc;
} Poly;

typedef struct {

  // Block statistics
  BlockSums blocks[NUM_BLOCKS];

  // Current block accumulators
  int64_t acc0;
  int64_t acc1;
  int64_t acc2;

  uint32_t u;
  uint32_t k;

  // Reading statistics
  uint32_t totalReadings;
  uint32_t worstCycles;

  // Glitch detection
  int32_t prevGood;
  uint32_t havePrev;
  uint32_t glitches;

  // Flight regions
  uint32_t coastBlocks;
  uint32_t lineFirstBlock;
  uint32_t ok;

  // Fitted curves
  Poly para;
  Poly line;

  // Precomputed powers
  double up[5];

  // Final results
  double gravity;
  double apogeeT;
  double apogeeH;
  double deployT;
  double deployH;
  double descent;
  double landT;

  // Serial command buffer
  char cmdBuf[40];
  uint8_t cmdLen;

} State;

static State state;

_Static_assert(sizeof(State) <= 4096,
               "State struct exceeds 4096 bytes!");

// ============================================================
// TRACKER
// ============================================================

static uint32_t trackerSeed = 20261001u;

static uint32_t trackerRand(void) {
  trackerSeed = trackerSeed * 1664525u + 1013904223u;
  return trackerSeed;
}

float nextAltitude(void) {

  const int64_t H0 = 31200;
  const int64_t V0 = 7850;
  const int64_t G  = 981;
  const int64_t TD = 12000;
  const int64_t VD = -450;

  static int64_t t = 0;
  static int quiet = 0;

  int64_t tc = (t < TD) ? t : TD;

  int64_t h =
      H0
      + V0 * tc / 1000
      - G * tc * tc / 2000000;

  if (t > TD) {
    h += VD * (t - TD) / 1000;
  }

  uint32_t r = trackerRand();

  if (quiet > 0) {
    quiet--;
  }

  if (TRACKER_GLITCHES &&
      t >= 1000 &&
      quiet == 0 &&
      (r >> 24) < 2) {

    quiet = 10;

    int64_t off =
        5000 + (int64_t)((r >> 4) % 25001);

    h += ((r >> 3) & 1) ? off : -off;

    if (h < 0) {
      h += 2 * off;
    }

  } else {

    h += (int64_t)((r >> 8) % 401) - 200;
  }

  t++;

  return (float)h / 100.0f;
}

// ============================================================
// ADD READING
// ============================================================

static inline void addReading(int32_t h) {

  int64_t hh = h;

  state.acc0 += hh;

  state.acc1 +=
      (int64_t)state.u * hh;

  state.acc2 +=
      (int64_t)(state.u * state.u) * hh;

  state.u++;

  if (state.u == BLOCK_LEN) {

    if (state.k < NUM_BLOCKS) {

      BlockSums *b =
          &state.blocks[state.k];

      b->s0 = state.acc0;
      b->s1 = state.acc1;
      b->s2 = state.acc2;
    }

    state.acc0 = 0;
    state.acc1 = 0;
    state.acc2 = 0;

    state.u = 0;
    state.k++;
  }
}

// ============================================================
// PRECOMPUTE LOCAL POWERS
// ============================================================

static void initLocalPowers(void) {

  double m =
      (double)BLOCK_LEN - 1.0;

  state.up[0] =
      (double)BLOCK_LEN;

  state.up[1] =
      m * (m + 1.0) / 2.0;

  state.up[2] =
      m * (m + 1.0) *
      (2.0 * m + 1.0) / 6.0;

  state.up[3] =
      state.up[1] * state.up[1];

  state.up[4] =
      m * (m + 1.0) *
      (2.0 * m + 1.0) *
      (3.0 * m * m +
       3.0 * m - 1.0) / 30.0;
}

// ============================================================
// POLYNOMIAL EVALUATION
// ============================================================

static double evalPoly(
    const Poly *p,
    double x) {

  double d = x - p->xc;

  return p->c0
       + p->c1 * d
       + p->c2 * d * d;
}

// ============================================================
// GAUSSIAN ELIMINATION
// ============================================================

static bool gaussSolve(
    int n,
    double M[3][4],
    double *sol) {

  for (int i = 0; i < n; i++) {

    int piv = i;

    for (int j = i + 1; j < n; j++) {

      if (fabs(M[j][i]) >
          fabs(M[piv][i])) {

        piv = j;
      }
    }

    if (fabs(M[piv][i]) < 1e-18) {
      return false;
    }

    if (piv != i) {

      for (int k = i; k <= n; k++) {

        double tmp = M[i][k];

        M[i][k] = M[piv][k];

        M[piv][k] = tmp;
      }
    }

    for (int j = i + 1; j < n; j++) {

      double f =
          M[j][i] / M[i][i];

      for (int k = i; k <= n; k++) {

        M[j][k] -=
            f * M[i][k];
      }
    }
  }

  for (int i = n - 1; i >= 0; i--) {

    double s = M[i][n];

    for (int j = i + 1; j < n; j++) {

      s -=
          M[i][j] * sol[j];
    }

    sol[i] =
        s / M[i][i];
  }

  return true;
}

// ============================================================
// FIT REGION
// ============================================================

static bool fitRegion(
    uint32_t kFrom,
    uint32_t kTo,
    int deg,
    Poly *out) {

  if (kTo <= kFrom) {
    return false;
  }

  double xcMs =
      ((double)kFrom * BLOCK_LEN
       + (double)kTo * BLOCK_LEN
       - 1.0) / 2.0;

  double sp[5] =
      {0, 0, 0, 0, 0};

  double sq[3] =
      {0, 0, 0};

  for (uint32_t k = kFrom;
       k < kTo;
       k++) {

    double o =
        (double)k * BLOCK_LEN
        - xcMs;

    double op[5];

    op[0] = 1.0;

    for (int i = 1; i < 5; i++) {
      op[i] =
          op[i - 1] * o;
    }

    double s[3];

    s[0] =
        (double)state.blocks[k].s0;

    s[1] =
        (double)state.blocks[k].s1;

    s[2] =
        (double)state.blocks[k].s2;

    for (int p = 0;
         p <= 2 * deg;
         p++) {

      double acc = 0.0;

      for (int j = 0;
           j <= p;
           j++) {

        acc +=
            BINOM[p][j]
            * op[p - j]
            * state.up[j];
      }

      sp[p] += acc;
    }

    for (int j = 0;
         j <= deg;
         j++) {

      double acc = 0.0;

      for (int i = 0;
           i <= j;
           i++) {

        acc +=
            BINOM[j][i]
            * op[j - i]
            * s[i];
      }

      sq[j] += acc;
    }
  }

  int n = deg + 1;

  double M[3][4];

  for (int i = 0;
       i < n;
       i++) {

    for (int j = 0;
         j < n;
         j++) {

      M[i][j] =
          sp[i + j]
          / K1000[i + j];
    }

    M[i][n] =
        sq[i]
        / (K1000[i] * 100.0);
  }

  double sol[3] =
      {0, 0, 0};

  if (!gaussSolve(
          n,
          M,
          sol)) {

    return false;
  }

  out->xc =
      xcMs / 1000.0;

  out->c0 =
      sol[0];

  out->c1 =
      sol[1];

  out->c2 =
      (deg == 2)
      ? sol[2]
      : 0.0;

  return true;
}

// ============================================================
// BLOCK SLOPE
// ============================================================

static double blockSlope(
    uint32_t k) {

  double n =
      (double)BLOCK_LEN;

  double num =
      n * (double)state.blocks[k].s1
      - state.up[1]
      * (double)state.blocks[k].s0;

  double den =
      n * state.up[2]
      - state.up[1]
      * state.up[1];

  return num / den * 10.0;
}

// ============================================================
// BLOCK MEAN
// ============================================================

static double blockMean(
    uint32_t k) {

  return
      (double)state.blocks[k].s0
      / (double)BLOCK_LEN
      / 100.0;
}

// ============================================================
// BLOCK MIDPOINT
// ============================================================

static double blockMidX(
    uint32_t k) {

  return
      ((double)k * BLOCK_LEN
       + ((double)BLOCK_LEN - 1.0)
         / 2.0)
      / 1000.0;
}

// ============================================================
// ANALYSE FLIGHT
// ============================================================

static void analyse(void) {

  state.ok = 0;

  initLocalPowers();

  Poly seed;

  if (!fitRegion(
          SURE_PARACHUTE_BLOCK,
          NUM_BLOCKS,
          1,
          &seed)) {

    return;
  }

  int kf =
      (int)SURE_PARACHUTE_BLOCK - 1;

  while (kf >= 0) {

    double dSlope =
        fabs(
            blockSlope(kf)
            - seed.c1);

    double dHeight =
        fabs(
            blockMean(kf)
            - evalPoly(
                &seed,
                blockMidX(kf)));

    if (dSlope > 1.2 ||
        dHeight > 0.4) {

      break;
    }

    kf--;
  }

  if (kf < 3) {
    return;
  }

  // kf is the transition block containing both coast and parachute data.
  // Exclude it from the coast fit so the parabola uses 100% pure coast data.
  state.coastBlocks =
      (kf > 3) ? (uint32_t)kf : 3;

  state.lineFirstBlock =
      (uint32_t)kf + 1;

  if (!fitRegion(
          state.lineFirstBlock,
          NUM_BLOCKS,
          1,
          &state.line)) {

    return;
  }

  if (!fitRegion(
          0,
          state.coastBlocks,
          2,
          &state.para)) {

    return;
  }

  const Poly *P =
      &state.para;

  const Poly *L =
      &state.line;

  // Descent speed
  state.descent =
      L->c1;

  // Landing time
  state.landT =
      L->xc
      - L->c0 / L->c1;

  // Gravity
  state.gravity =
      -2.0 * P->c2;

  // Apogee
  state.apogeeT =
      P->xc
      - P->c1
        / (2.0 * P->c2);

  state.apogeeH =
      P->c0
      - P->c1 * P->c1
        / (4.0 * P->c2);

  // Parabola - line intersection
  double qa =
      P->c2;

  double qb =
      P->c1
      - 2.0 * P->c2
        * P->xc;

  double qc =
      P->c0
      - P->c1 * P->xc
      + P->c2
        * P->xc
        * P->xc
      - (L->c0
         - L->c1
           * L->xc);

  double q_lin_a =
      qa;

  double q_lin_b =
      qb - L->c1;

  double q_lin_c =
      qc;

  double disc =
      q_lin_b * q_lin_b
      - 4.0
        * q_lin_a
        * q_lin_c;

  if (disc < 0.0) {
    return;
  }

  double root =
      sqrt(disc);

  double r1 =
      (-q_lin_b + root)
      / (2.0 * q_lin_a);

  double r2 =
      (-q_lin_b - root)
      / (2.0 * q_lin_a);

  state.deployT =
      (r1 > r2)
      ? r1
      : r2;

  state.deployH =
      evalPoly(
          P,
          state.deployT);

  state.ok = 1;
}

// ============================================================
// ALTITUDE AT TIME
// ============================================================

static double altitudeAt(
    double t) {

  if (t < state.deployT) {

    return evalPoly(
        &state.para,
        t);

  }

  return evalPoly(
      &state.line,
      t);
}

// ============================================================
// PROCESS ALL READINGS
// ============================================================

static void processReadings(void) {

  for (uint32_t i = 0;
       i < TOTAL_READINGS;
       i++) {

    // Tracker execution is intentionally
    // outside the timing measurement.
    float altitude =
        nextAltitude();

    uint32_t start =
        ESP.getCycleCount();

    int32_t h =
        (int32_t)lroundf(
            altitude * 100.0f);

    // Glitch filtering
    if (state.havePrev) {

      int32_t d =
          h - state.prevGood;

      if (d < 0) {
        d = -d;
      }

      if (d > GLITCH_JUMP_CM) {

        state.glitches++;

        h =
            state.prevGood;

      } else {

        state.prevGood =
            h;
      }

    } else {

      state.prevGood =
          h;

      state.havePrev =
          1;
    }

    addReading(h);

    state.totalReadings++;

    uint32_t cycles =
        ESP.getCycleCount()
        - start;

    if (cycles >
        state.worstCycles) {

      state.worstCycles =
          cycles;
    }
  }

  analyse();
}

// ============================================================
// PRINT COMMAND PROMPT
// ============================================================

static void printPrompt(void) {

  Serial.println();
  Serial.println("----------------------------------------");
  Serial.println("Enter command:");
  Serial.print("> ");
}

// ============================================================
// HANDLE COMMAND
// ============================================================

static void handleCommand(
    char *cmd) {

  while (*cmd == ' ') {
    cmd++;
  }

  int len =
      (int)strlen(cmd);

  while (len > 0 &&
         (cmd[len - 1] == ' ' ||
          cmd[len - 1] == '\r')) {

    cmd[--len] = 0;
  }

  // Convert command to uppercase
  for (int i = 0;
       i < len;
       i++) {

    if (cmd[i] >= 'a' &&
        cmd[i] <= 'z') {

      cmd[i] -= 32;
    }
  }

  uint32_t start =
      ESP.getCycleCount();

  // Commands requiring valid analysis
  if (strcmp(cmd, "STATS") != 0 &&
      strcmp(cmd, "GLITCHES") != 0 &&
      !state.ok) {

    Serial.println("ERR");
    printPrompt();
    return;
  }

  // ----------------------------------------------------------
  // DESCENT
  // ----------------------------------------------------------

  if (strcmp(cmd, "DESCENT") == 0) {

    uint32_t cyc =
        ESP.getCycleCount()
        - start;

    Serial.print("DESCENT = ");
    Serial.print(
        state.descent,
        4);

    Serial.print(" m/s");
    Serial.print("    cycles = ");
    Serial.println(cyc);

    printPrompt();
    return;
  }

  // ----------------------------------------------------------
  // LAND
  // ----------------------------------------------------------

  if (strcmp(cmd, "LAND") == 0) {

    uint32_t cyc =
        ESP.getCycleCount()
        - start;

    Serial.print("LAND = ");
    Serial.print(
        state.landT,
        4);

    Serial.print(" s");
    Serial.print("    cycles = ");
    Serial.println(cyc);

    printPrompt();
    return;
  }

  // ----------------------------------------------------------
  // GRAVITY
  // ----------------------------------------------------------

  if (strcmp(cmd, "G") == 0) {

    uint32_t cyc =
        ESP.getCycleCount()
        - start;

    Serial.print("G = ");
    Serial.print(
        state.gravity,
        4);

    Serial.print(" m/s^2");
    Serial.print("    cycles = ");
    Serial.println(cyc);

    printPrompt();
    return;
  }

  // ----------------------------------------------------------
  // APOGEE
  // ----------------------------------------------------------

  if (strcmp(cmd, "APOGEE") == 0) {

    uint32_t cyc =
        ESP.getCycleCount()
        - start;

    Serial.print("APOGEE = ");
    Serial.print(
        state.apogeeT,
        4);

    Serial.print(" s, height = ");

    Serial.print(
        state.apogeeH,
        4);

    Serial.print(" m");
    Serial.print("    cycles = ");
    Serial.println(cyc);

    printPrompt();
    return;
  }

  // ----------------------------------------------------------
  // DEPLOY
  // ----------------------------------------------------------

  if (strcmp(cmd, "DEPLOY") == 0) {

    uint32_t cyc =
        ESP.getCycleCount()
        - start;

    Serial.print("DEPLOY = ");
    Serial.print(
        state.deployT,
        4);

    Serial.print(" s, altitude = ");

    Serial.print(
        state.deployH,
        4);

    Serial.print(" m");
    Serial.print("    cycles = ");
    Serial.println(cyc);

    printPrompt();
    return;
  }

  // ----------------------------------------------------------
  // GLITCHES
  // ----------------------------------------------------------

  if (strcmp(cmd, "GLITCHES") == 0) {

    uint32_t cyc =
        ESP.getCycleCount()
        - start;

    Serial.print("GLITCHES = ");
    Serial.print(
        state.glitches);

    Serial.print("    cycles = ");
    Serial.println(cyc);

    printPrompt();
    return;
  }

  // ----------------------------------------------------------
  // ALTITUDE
  // ----------------------------------------------------------

  if (strncmp(cmd, "ALT ", 4) == 0) {

    char *endPtr;

    double t =
        strtod(
            cmd + 4,
            &endPtr);

    if (endPtr == cmd + 4 ||
        *endPtr != '\0' ||
        t < 0.0 ||
        t > 100.0) {

      Serial.println("ERR");
      printPrompt();
      return;
    }

    double alt =
        altitudeAt(t);

    uint32_t cyc =
        ESP.getCycleCount()
        - start;

    Serial.print("ALT ");
    Serial.print(t, 3);

    Serial.print(" s = ");

    Serial.print(
        alt,
        3);

    Serial.print(" m");
    Serial.print("    cycles = ");
    Serial.println(cyc);

    printPrompt();
    return;
  }

  // ----------------------------------------------------------
  // STATS
  // ----------------------------------------------------------

  if (strcmp(cmd, "STATS") == 0) {

    uint32_t cyc =
        ESP.getCycleCount()
        - start;

    Serial.println();
    Serial.println("SYSTEM STATS");
    Serial.println("----------------------------------------");

    Serial.print("Memory       : ");
    Serial.print(sizeof(State));
    Serial.println(" bytes / 4096");

    Serial.print("Worst cycles : ");
    Serial.print(state.worstCycles);
    Serial.println(" / 12000");

    Serial.print("Readings     : ");
    Serial.println(state.totalReadings);

    Serial.print("Glitches     : ");
    Serial.println(state.glitches);

    Serial.print("Coast blocks : ");
    Serial.println(state.coastBlocks);

    Serial.print("Line block   : ");
    Serial.println(state.lineFirstBlock);

    Serial.print("Command cycles: ");
    Serial.println(cyc);

    printPrompt();
    return;
  }

  // ----------------------------------------------------------
  // UNKNOWN COMMAND
  // ----------------------------------------------------------

  Serial.println("ERR");
  printPrompt();
}

// ============================================================
// SERIAL SERVICE
// ============================================================

static void serviceSerial(void) {

  while (Serial.available()) {

    char c =
        (char)Serial.read();

    if (c == '\r' ||
        c == '\n') {

      if (state.cmdLen > 0) {

        state.cmdBuf[
            state.cmdLen] =
            '\0';

        handleCommand(
            state.cmdBuf);

        state.cmdLen = 0;
      }

      continue;
    }

    if (state.cmdLen <
        sizeof(state.cmdBuf) - 1) {

      state.cmdBuf[
          state.cmdLen++] =
          c;

    } else {

      state.cmdLen = 0;

      Serial.println("ERR");

      printPrompt();
    }
  }
}

// ============================================================
// SETUP
// ============================================================

void setup(void) {

  Serial.begin(115200);

  delay(1000);

  processReadings();

  Serial.println();
  Serial.println("========================================");
  Serial.println("           CATCH THE ROCKET");
  Serial.println("        Flight Analysis System");
  Serial.println("========================================");

  Serial.println();
  Serial.println("READY");

  Serial.println();
  Serial.print("MEMORY       : ");
  Serial.print(sizeof(State));
  Serial.println(" bytes / 4096");

  Serial.print("WORST CYCLES : ");
  Serial.print(state.worstCycles);
  Serial.println(" / 12000");

  Serial.println();
  Serial.println("Available commands:");
  Serial.println("  G");
  Serial.println("  APOGEE");
  Serial.println("  DEPLOY");
  Serial.println("  DESCENT");
  Serial.println("  LAND");
  Serial.println("  GLITCHES");
  Serial.println("  ALT <time>");
  Serial.println("  STATS");

  printPrompt();
}

// ============================================================
// LOOP
// ============================================================

void loop(void) {

  serviceSerial();
}