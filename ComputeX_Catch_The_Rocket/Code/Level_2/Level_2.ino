#define TRACKER_GLITCHES 0

// ===== OFFICIAL TRACKER: do not modify =====
static uint32_t trackerSeed = 20261001u;

static uint32_t trackerRand() {
  trackerSeed = trackerSeed * 1664525u + 1013904223u;
  return trackerSeed;
}

float nextAltitude() {
  const int64_t H0 = 31200;
  const int64_t V0 = 7850;
  const int64_t G = 981;
  const int64_t TD = 12000;
  const int64_t VD = -450;

  static int64_t t = 0;
  static int quiet = 0;

  int64_t tc = (t < TD) ? t : TD;

  int64_t h = H0
            + V0 * tc / 1000
            - G * tc * tc / 2000000;

  if (t > TD)
    h += VD * (t - TD) / 1000;

  uint32_t r = trackerRand();

  if (TRACKER_GLITCHES && t >= 1000 && quiet == 0 && (r >> 24) < 2) {
    quiet = 10;

    int64_t off = 5000 + (int64_t)((r >> 4) % 25001);

    h += ((r >> 3) & 1) ? off : -off;

    if (h < 0)
      h += 2 * off;
  }
  else {
    h += (int64_t)((r >> 8) % 401) - 200;
  }

  t++;

  return (float)h / 100.0f;
}
// ===== END OFFICIAL TRACKER =====


// ===== STATE =====

struct State {
  uint32_t n;

  int64_t sumT;
  int64_t sumT2;
  int64_t sumT3;

  uint64_t sumT4;

  int64_t sumH;
  int64_t sumTH;
  int64_t sumT2H;

  uint32_t worstCycles;

  double gravity;
  double apogeeTime;
  double apogeeHeight;
};

State state = {0};


// ===== SOLVE 3 x 3 SYSTEM =====
//
// [ n     sumT   sumT2 ] [a]   [sumH]
// [ sumT  sumT2  sumT3 ] [b] = [sumTH]
// [ sumT2 sumT3  sumT4 ] [c]   [sumT2H]
//
// Gaussian elimination is used.

void solveQuadratic() {

  double A[3][4];

  A[0][0] = (double)state.n;
  A[0][1] = (double)state.sumT;
  A[0][2] = (double)state.sumT2;
  A[0][3] = (double)state.sumH;

  A[1][0] = (double)state.sumT;
  A[1][1] = (double)state.sumT2;
  A[1][2] = (double)state.sumT3;
  A[1][3] = (double)state.sumTH;

  A[2][0] = (double)state.sumT2;
  A[2][1] = (double)state.sumT3;
  A[2][2] = (double)state.sumT4;
  A[2][3] = (double)state.sumT2H;


  // Forward elimination

  for (int i = 0; i < 3; i++) {

    double pivot = A[i][i];

    for (int j = i + 1; j < 3; j++) {

      double factor = A[j][i] / pivot;

      for (int k = i; k < 4; k++) {
        A[j][k] -= factor * A[i][k];
      }
    }
  }


  // Back substitution

  double c = A[2][3] / A[2][2];

  double b =
      (A[1][3] - A[1][2] * c) / A[1][1];

  double a =
      (A[0][3]
       - A[0][2] * c
       - A[0][1] * b) / A[0][0];


  // h = a + b*t + c*t²
  //
  // t is in milliseconds.
  // h is in centimetres.

  // g:
  // c = -g / 2, after unit conversion.
  // Therefore:
  // g(m/s²) = -2*c*10000

  state.gravity = -2.0 * c * 10000.0;


  // Apogee occurs where dh/dt = 0:
  //
  // b + 2ct = 0
  //
  // t = -b/(2c)

  double tApex = -b / (2.0 * c);

  state.apogeeTime = tApex / 1000.0;

  state.apogeeHeight =
      (a + b * tApex + c * tApex * tApex) / 100.0;
}


// ===== PROCESS 100,000 READINGS =====

void processReadings() {

  for (uint32_t i = 0; i < 100000; i++) {

    uint32_t start = ESP.getCycleCount();

    float altitude = nextAltitude();


    // First 8 seconds = 8000 readings
    if (i < 8000) {

      // milliseconds
      int64_t t = (int64_t)i;

      // metres -> centimetres
      int64_t h =
          (int64_t)(altitude * 100.0);


      state.n++;

      state.sumT += t;

      state.sumT2 += t * t;

      state.sumT3 += t * t * t;

      state.sumT4 +=
          (uint64_t)t * t * t * t;

      state.sumH += h;

      state.sumTH += t * h;

      state.sumT2H += t * t * h;
    }


    uint32_t end = ESP.getCycleCount();

    uint32_t cycles = end - start;

    if (cycles > state.worstCycles)
      state.worstCycles = cycles;
  }


  solveQuadratic();
}


// ===== COMMAND HANDLER =====

void handleCommand(String command) {

  command.trim();


  if (command == "G") {

    uint32_t start = ESP.getCycleCount();

    double answer = state.gravity;

    uint32_t cycles =
        ESP.getCycleCount() - start;

    Serial.print("G = ");
    Serial.print(answer, 4);
    Serial.print(" m/s^2, cycles = ");
    Serial.println(cycles);
  }


  else if (command == "APOGEE") {

    uint32_t start = ESP.getCycleCount();

    double time = state.apogeeTime;
    double height = state.apogeeHeight;

    uint32_t cycles =
        ESP.getCycleCount() - start;

    Serial.print("APOGEE time = ");
    Serial.print(time, 4);
    Serial.print(" s, height = ");
    Serial.print(height, 4);
    Serial.print(" m, cycles = ");
    Serial.println(cycles);
  }


  else if (command == "STATS") {

    uint32_t start = ESP.getCycleCount();

    Serial.print("Readings used = ");
    Serial.println(state.n);

    Serial.print("State size = ");
    Serial.print(sizeof(State));
    Serial.println(" bytes");

    Serial.print("Worst cycles = ");
    Serial.println(state.worstCycles);

    uint32_t cycles =
        ESP.getCycleCount() - start;

    Serial.print("STATS cycles = ");
    Serial.println(cycles);
  }


  else {

    Serial.println("ERR");
  }
}


// ===== SETUP =====

void setup() {

  Serial.begin(115200);

  Serial.println("Catch the Rocket - Level 2");

  Serial.print("State size = ");
  Serial.print(sizeof(State));
  Serial.println(" bytes");

  processReadings();

  Serial.println("Ready for commands.");
}


// ===== LOOP =====

void loop() {

  if (Serial.available()) {

    String command =
        Serial.readStringUntil('\n');

    handleCommand(command);
  }
}