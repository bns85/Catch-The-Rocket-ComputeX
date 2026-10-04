// ============================================================
// Catch the Rocket - Level 3 (block-sum design)
//
// IDEA
//   * The 100,000 readings are cut into 100 blocks of 1,000 readings.
//   * For each block we keep only THREE exact integers (using the
//     reading's position u = 0..999 inside the block, h in cm):
//         s0 = sum(h)     s1 = sum(u*h)     s2 = sum(u*u*h)
//     The time grid is the same in every block, so the sums that only
//     involve time (sum u, sum u^2, ...) are known formulas, not stored.
//   * After the burst, sums over any group of blocks are rebuilt in
//     double precision, so we can fit a line or a parabola to ANY
//     group of blocks without having kept a single reading.
//   * The parachute is found from the DATA (no 12 s hard-coding):
//       1. fit a line to everything from 45 s on (the problem says the
//          rocket is under the parachute by then),
//       2. walk backwards block by block while a block still agrees
//          with that line (slope AND height),
//       3. the first block that disagrees is dropped; everything
//          before it is the coast -> fit the parabola to it.
//   * DEPLOY = where parabola and line cross. The parabola minus the
//     line opens downward, so it has two roots. At the LARGER root the
//     parabola goes from above the line to below it (it falls faster
//     than the parachute line, by at least 10 m/s). That is the real
//     deployment; the smaller root is during the climb.
//
// MEMORY  : 100 blocks x 24 bytes + a few results = well under 4096.
// SPEED   : per reading: 1 float->int, 3 adds, 2 multiplies.
// ============================================================

#include <Arduino.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define TRACKER_GLITCHES 0

// NOTE: the types and constants come first so the Arduino IDE's automatic
// function prototypes (which it inserts before the first function) can see them.

// ============================================================
// DESIGN CONSTANTS
// ============================================================

static const uint32_t TOTAL_READINGS = 100000;
static const uint32_t BLOCK_LEN      = 1000;                   // readings per block
static const uint32_t NUM_BLOCKS     = TOTAL_READINGS / BLOCK_LEN;

// The problem guarantees the parachute is open from 45 s onward (Level 1).
static const uint32_t SURE_PARACHUTE_BLOCK = 45000 / BLOCK_LEN;

// A block "agrees with the parachute line" if its slope is within
// 1 m/s of the line (block-slope noise is only about 0.13 m/s) AND its
// mean height is within 0.3 m of the line (noise about 0.04 m). A coast
// block next to the deployment fails at least one of the two tests, even
// for the nastiest flight the 10 m/s guarantee allows.
static const double SLOPE_TOL_MPS = 1.0;
static const double HEIGHT_TOL_M  = 0.3;

// A glitch moves one reading 50 m to 300 m, while two normal readings
// differ by at most 4 m plus the motion in 1-2 ms. Anything more than
// 25 m away from the last good reading is a glitch. The rocket's
// altitude is continuous when the parachute opens (only the SLOPE
// changes), so a glitch can never be mistaken for the deployment.
static const int32_t GLITCH_JUMP_CM = 2500;

// binomial coefficients and powers of 1000 (read-only, live in flash)
static const double BINOM[5][5] = {
  {1, 0, 0, 0, 0},
  {1, 1, 0, 0, 0},
  {1, 2, 1, 0, 0},
  {1, 3, 3, 1, 0},
  {1, 4, 6, 4, 1}
};
static const double K1000[5] = {1.0, 1e3, 1e6, 1e9, 1e12};


// ============================================================
// STATE (everything the board remembers)
// ============================================================

struct BlockSums {
  int64_t s0;   // sum of h          (h in cm)
  int64_t s1;   // sum of u * h      (u = position in block, ms)
  int64_t s2;   // sum of u * u * h
};

// y = c0 + c1*(x - xc) + c2*(x - xc)^2     x in seconds, y in metres
struct Poly {
  double c0, c1, c2, xc;
};

struct State {
  BlockSums blocks[NUM_BLOCKS];

  // the block currently being filled
  int64_t  acc0, acc1, acc2;
  uint32_t u;                 // position inside current block
  uint32_t k;                 // index of current block

  uint32_t totalReadings;
  uint32_t worstCycles;

  // glitch filter
  int32_t  prevGood;          // last reading accepted as good (cm)
  uint32_t havePrev;
  uint32_t glitches;          // exact count of rejected readings

  // what the data told us
  uint32_t coastBlocks;       // blocks 0 .. coastBlocks-1 are pure coast
  uint32_t lineFirstBlock;    // blocks lineFirstBlock .. 99 are on the parachute line
  uint32_t ok;                // 1 if the analysis succeeded

  Poly para;                  // coast parabola
  Poly line;                  // parachute line (c2 = 0)

  double up[5];               // sum of u^j over one block, j = 0..4

  // answers
  double gravity;             // m/s^2
  double apogeeT, apogeeH;    // s, m
  double deployT, deployH;    // s, m
  double descent;             // m/s (negative = falling)
  double landT;               // s
};

State state;

// The compiler refuses to build the sketch if the remembered data ever exceeds 4 KB.
static_assert(sizeof(State) <= 4096, "State is larger than 4096 bytes");



// ===== OFFICIAL TRACKER: do not modify =====
// Set TRACKER_GLITCHES to 1 or 0 *above* this block, never inside it.
static uint32_t trackerSeed = 20261001u;

static uint32_t trackerRand() {
  trackerSeed = trackerSeed * 1664525u + 1013904223u;
  return trackerSeed;
}

float nextAltitude() {
  // This flight's secrets. Every tracker uses different ones.
  const int64_t H0 = 31200;   // altitude at t = 0, in cm
  const int64_t V0 = 7850;    // upward speed at t = 0, in cm/s
  const int64_t G = 981;      // gravity, in cm/s^2
  const int64_t TD = 12000;   // parachute opens at this time, in ms
  const int64_t VD = -450;    // speed under the parachute, in cm/s
  static int64_t t = 0;       // time of this reading, in ms
  static int quiet = 0;       // readings until a glitch is allowed again
  int64_t tc = (t < TD) ? t : TD;
  int64_t h = H0 + V0 * tc / 1000 - G * tc * tc / 2000000;
  if (t > TD) h += VD * (t - TD) / 1000;
  uint32_t r = trackerRand();
  if (quiet > 0) quiet--;
  if (TRACKER_GLITCHES && t >= 1000 && quiet == 0 && (r >> 24) < 2) {
    quiet = 10;
    // glitch: 50 to 300 m off
    int64_t off = 5000 + (int64_t)((r >> 4) % 25001);
    h += ((r >> 3) & 1) ? off : -off;
    if (h < 0) h += 2 * off;
  } else {
    h += (int64_t)((r >> 8) % 401) - 200;   // noise: -2.00 to +2.00 m
  }
  t++;
  return (float)h / 100.0f;
}
// ========================================= END OFFICIAL TRACKER

// ============================================================
// PER-READING WORK (the only code that runs inside the burst)
// ============================================================

static inline void addReading(int32_t h) {
  int64_t hh = h;
  state.acc0 += hh;
  state.acc1 += (int64_t)state.u * hh;
  state.acc2 += (int64_t)(state.u * state.u) * hh;

  if (++state.u == BLOCK_LEN) {
    if (state.k < NUM_BLOCKS) {
      BlockSums &b = state.blocks[state.k];
      b.s0 = state.acc0;
      b.s1 = state.acc1;
      b.s2 = state.acc2;
    }
    state.acc0 = state.acc1 = state.acc2 = 0;
    state.u = 0;
    state.k++;
  }
}


// ============================================================
// MATHS AFTER THE BURST
// ============================================================

static void initLocalPowers() {
  double m = (double)BLOCK_LEN - 1.0;
  state.up[0] = (double)BLOCK_LEN;
  state.up[1] = m * (m + 1.0) / 2.0;
  state.up[2] = m * (m + 1.0) * (2.0 * m + 1.0) / 6.0;
  state.up[3] = state.up[1] * state.up[1];
  state.up[4] = m * (m + 1.0) * (2.0 * m + 1.0) * (3.0 * m * m + 3.0 * m - 1.0) / 30.0;
}

static double evalPoly(const Poly &p, double x) {
  double d = x - p.xc;
  return p.c0 + p.c1 * d + p.c2 * d * d;
}

// Solve an n x n system (n <= 3) by Gaussian elimination with pivoting.
static bool gaussSolve(int n, double M[3][4], double *sol) {
  for (int i = 0; i < n; i++) {
    int piv = i;
    for (int j = i + 1; j < n; j++)
      if (fabs(M[j][i]) > fabs(M[piv][i])) piv = j;
    if (fabs(M[piv][i]) < 1e-18) return false;
    if (piv != i)
      for (int k = i; k <= n; k++) {
        double tmp = M[i][k]; M[i][k] = M[piv][k]; M[piv][k] = tmp;
      }
    for (int j = i + 1; j < n; j++) {
      double f = M[j][i] / M[i][i];
      for (int k = i; k <= n; k++) M[j][k] -= f * M[i][k];
    }
  }
  for (int i = n - 1; i >= 0; i--) {
    double s = M[i][n];
    for (int j = i + 1; j < n; j++) s -= M[i][j] * sol[j];
    sol[i] = s / M[i][i];
  }
  return true;
}

// Least-squares polynomial (degree 1 or 2) over blocks kFrom .. kTo-1.
// Time is measured from the middle of the region (in seconds) so all
// numbers stay small and well-conditioned.
static bool fitRegion(uint32_t kFrom, uint32_t kTo, int deg, Poly &out) {
  double xcMs = ((double)kFrom * BLOCK_LEN + (double)kTo * BLOCK_LEN - 1.0) / 2.0;

  double sp[5] = {0, 0, 0, 0, 0};   // sum of x^p   (ms for now)
  double sq[3] = {0, 0, 0};         // sum of x^j * h

  for (uint32_t k = kFrom; k < kTo; k++) {
    double o = (double)k * BLOCK_LEN - xcMs;      // block start relative to centre
    double op[5];
    op[0] = 1.0;
    for (int i = 1; i < 5; i++) op[i] = op[i - 1] * o;

    double s[3] = {(double)state.blocks[k].s0,
                   (double)state.blocks[k].s1,
                   (double)state.blocks[k].s2};

    // sum over the block of (u + o)^p
    for (int p = 0; p <= 2 * deg; p++) {
      double acc = 0;
      for (int j = 0; j <= p; j++) acc += BINOM[p][j] * op[p - j] * state.up[j];
      sp[p] += acc;
    }
    // sum over the block of (u + o)^j * h
    for (int j = 0; j <= deg; j++) {
      double acc = 0;
      for (int i = 0; i <= j; i++) acc += BINOM[j][i] * op[j - i] * s[i];
      sq[j] += acc;
    }
  }

  int n = deg + 1;
  double M[3][4];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) M[i][j] = sp[i + j] / K1000[i + j];       // seconds
    M[i][n] = sq[i] / (K1000[i] * 100.0);                                 // cm -> m
  }

  double sol[3] = {0, 0, 0};
  if (!gaussSolve(n, M, sol)) return false;

  out.xc = xcMs / 1000.0;
  out.c0 = sol[0];
  out.c1 = sol[1];
  out.c2 = (deg == 2) ? sol[2] : 0.0;
  return true;
}

// slope (m/s) of the straight line fitted to one block
static double blockSlope(uint32_t k) {
  double n   = (double)BLOCK_LEN;
  double num = n * (double)state.blocks[k].s1 - state.up[1] * (double)state.blocks[k].s0;
  double den = n * state.up[2] - state.up[1] * state.up[1];
  return num / den * 10.0;                       // cm/ms -> m/s
}

static double blockMean(uint32_t k) {            // metres
  return (double)state.blocks[k].s0 / (double)BLOCK_LEN / 100.0;
}

static double blockMidX(uint32_t k) {            // seconds
  return ((double)k * BLOCK_LEN + ((double)BLOCK_LEN - 1.0) / 2.0) / 1000.0;
}

static void analyse() {
  state.ok = 0;
  initLocalPowers();

  // 1. seed line from the part we are told is under the parachute
  Poly seed;
  if (!fitRegion(SURE_PARACHUTE_BLOCK, NUM_BLOCKS, 1, seed)) return;

  // 2. walk backwards while blocks still agree with that line
  int kf = (int)SURE_PARACHUTE_BLOCK - 1;
  while (kf >= 0) {
    double dSlope  = fabs(blockSlope(kf) - seed.c1);
    double dHeight = fabs(blockMean(kf) - evalPoly(seed, blockMidX(kf)));
    if (dSlope > SLOPE_TOL_MPS || dHeight > HEIGHT_TOL_M) break;
    kf--;
  }
  if (kf < 3) return;                            // not enough coast to fit a parabola

  // block kf is dropped: it is either pure coast or contains the deployment
  state.coastBlocks    = (uint32_t)kf;
  state.lineFirstBlock = (uint32_t)kf + 1;

  // 3. final fits on clean data
  if (!fitRegion(state.lineFirstBlock, NUM_BLOCKS, 1, state.line)) return;
  if (!fitRegion(0, state.coastBlocks, 2, state.para)) return;

  // 4. answers
  const Poly &P = state.para;
  const Poly &L = state.line;

  state.descent = L.c1;
  state.landT   = L.xc - L.c0 / L.c1;

  state.gravity = -2.0 * P.c2;
  state.apogeeT = P.xc - P.c1 / (2.0 * P.c2);
  state.apogeeH = P.c0 - P.c1 * P.c1 / (4.0 * P.c2);

  // parabola and line in absolute time, then solve  parabola = line
  double A2 = P.c2;
  double A1 = P.c1 - 2.0 * P.c2 * P.xc;
  double A0 = P.c0 - P.c1 * P.xc + P.c2 * P.xc * P.xc;
  double L1 = L.c1;
  double L0 = L.c0 - L.c1 * L.xc;

  double qa = A2, qb = A1 - L1, qc = A0 - L0;
  double disc = qb * qb - 4.0 * qa * qc;
  if (disc < 0.0) return;

  double sq = sqrt(disc);
  double r1 = (-qb + sq) / (2.0 * qa);
  double r2 = (-qb - sq) / (2.0 * qa);
  state.deployT = (r1 > r2) ? r1 : r2;           // the LARGER root (see header)
  state.deployH = evalPoly(P, state.deployT);

  state.ok = 1;
}

static double altitudeAt(double t) {
  return (t < state.deployT) ? evalPoly(state.para, t) : evalPoly(state.line, t);
}


// ============================================================
// RECEIVE THE BURST
// ============================================================

static void processReadings() {
  for (uint32_t i = 0; i < TOTAL_READINGS; i++) {

    uint32_t start = ESP.getCycleCount();

    float altitude = nextAltitude();
    int32_t h = (int32_t)lroundf(altitude * 100.0f);   // exact centimetres

    if (state.havePrev) {
      int32_t d = h - state.prevGood;
      if (d < 0) d = -d;
      if (d > GLITCH_JUMP_CM) {
        state.glitches++;
        h = state.prevGood;       // keep the time grid: reuse the last good value
      } else {
        state.prevGood = h;
      }
    } else {
      state.prevGood = h;
      state.havePrev = 1;
    }

    addReading(h);
    state.totalReadings++;

    uint32_t cycles = ESP.getCycleCount() - start;
    if (cycles > state.worstCycles) state.worstCycles = cycles;
  }

  analyse();
}


// ============================================================
// COMMANDS
// ============================================================

static void handleCommand(char *cmd) {

  // trim and upper-case in place
  while (*cmd == ' ') cmd++;
  int len = (int)strlen(cmd);
  while (len > 0 && (cmd[len - 1] == ' ' || cmd[len - 1] == '\r')) cmd[--len] = 0;
  for (int i = 0; i < len; i++)
    if (cmd[i] >= 'a' && cmd[i] <= 'z') cmd[i] -= 32;

  uint32_t start = ESP.getCycleCount();

  if (!state.ok && strcmp(cmd, "STATS") != 0 && strcmp(cmd, "GLITCHES") != 0) {
    Serial.println("ERR analysis failed");
    return;
  }

  if (strcmp(cmd, "DESCENT") == 0) {
    uint32_t cyc = ESP.getCycleCount() - start;
    Serial.print("DESCENT = ");
    Serial.print(state.descent, 4);
    Serial.print(" m/s, cycles = ");
    Serial.println(cyc);
  }
  else if (strcmp(cmd, "LAND") == 0) {
    uint32_t cyc = ESP.getCycleCount() - start;
    Serial.print("LAND = ");
    Serial.print(state.landT, 4);
    Serial.print(" s, cycles = ");
    Serial.println(cyc);
  }
  else if (strcmp(cmd, "G") == 0) {
    uint32_t cyc = ESP.getCycleCount() - start;
    Serial.print("G = ");
    Serial.print(state.gravity, 4);
    Serial.print(" m/s^2, cycles = ");
    Serial.println(cyc);
  }
  else if (strcmp(cmd, "APOGEE") == 0) {
    uint32_t cyc = ESP.getCycleCount() - start;
    Serial.print("APOGEE time = ");
    Serial.print(state.apogeeT, 4);
    Serial.print(" s, height = ");
    Serial.print(state.apogeeH, 4);
    Serial.print(" m, cycles = ");
    Serial.println(cyc);
  }
  else if (strcmp(cmd, "DEPLOY") == 0) {
    uint32_t cyc = ESP.getCycleCount() - start;
    Serial.print("DEPLOY time = ");
    Serial.print(state.deployT, 4);
    Serial.print(" s, altitude = ");
    Serial.print(state.deployH, 4);
    Serial.print(" m, cycles = ");
    Serial.println(cyc);
  }
  else if (strcmp(cmd, "GLITCHES") == 0) {
    uint32_t cyc = ESP.getCycleCount() - start;
    Serial.print("GLITCHES = ");
    Serial.print(state.glitches);
    Serial.print(", cycles = ");
    Serial.println(cyc);
  }
  else if (strncmp(cmd, "ALT ", 4) == 0) {
    char *endPtr;
    double t = strtod(cmd + 4, &endPtr);
    while (*endPtr == ' ') endPtr++;
    if (endPtr == cmd + 4 || *endPtr != '\0' || t < 0.0 || t > 100.0) {
      Serial.println("ERR");
      return;
    }
    double alt = altitudeAt(t);
    uint32_t cyc = ESP.getCycleCount() - start;
    Serial.print("ALT ");
    Serial.print(t, 3);
    Serial.print(" = ");
    Serial.print(alt, 2);
    Serial.print(" m, cycles = ");
    Serial.println(cyc);
  }
  else if (strcmp(cmd, "STATS") == 0) {
    Serial.print("Total readings = ");
    Serial.println(state.totalReadings);
    Serial.print("Coast readings used = ");
    Serial.println(state.coastBlocks * BLOCK_LEN);
    Serial.print("Parachute readings used = ");
    Serial.println((NUM_BLOCKS - state.lineFirstBlock) * BLOCK_LEN);
    Serial.print("Glitches rejected = ");
    Serial.println(state.glitches);
    Serial.print("Block length = ");
    Serial.println(BLOCK_LEN);
    Serial.print("State size = ");
    Serial.print(sizeof(State));
    Serial.println(" bytes");
    Serial.print("Worst cycles = ");
    Serial.println(state.worstCycles);
  }
  else {
    Serial.println("ERR");
  }
}


// ============================================================
// SETUP / LOOP
// ============================================================

void setup() {
  Serial.begin(115200);
  delay(1000);

  Serial.println();
  Serial.println("Catch the Rocket - Level 3");

  processReadings();

  Serial.println("Ready for commands.");
  Serial.print("State size = ");
  Serial.print(sizeof(State));
  Serial.println(" bytes");
}

void loop() {
  static char buf[32];
  static uint8_t n = 0;
  static uint32_t lastChar = 0;

  while (Serial.available()) {
    char ch = (char)Serial.read();
    lastChar = millis();
    if (ch == '\n' || ch == '\r') {
      if (n > 0) {
        buf[n] = 0;
        n = 0;
        handleCommand(buf);
      }
    } else if (n < sizeof(buf) - 1) {
      buf[n++] = ch;
    }
  }

  // Serial Monitor set to "No line ending": run the command after a short pause
  if (n > 0 && millis() - lastChar > 300) {
    buf[n] = 0;
    n = 0;
    handleCommand(buf);
  }
}
