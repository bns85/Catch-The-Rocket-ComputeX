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


// ===== OUR STATE =====

struct State {
  uint32_t n;

  int64_t sumT;
  int64_t sumH;
  int64_t sumTH;
  int64_t sumT2;

  uint32_t worstCycles;

  double descent;
  double landTime;
};

State state = {0, 0, 0, 0, 0, 0, 0, 0};


// ===== PROCESS 100,000 READINGS =====

void processReadings() {

  for (uint32_t i = 0; i < 100000; i++) {

    uint32_t start = ESP.getCycleCount();

    float altitude = nextAltitude();

    // Use readings from 45 seconds onward
    if (i >= 45000) {

      // 45 seconds becomes t' = 0
      int64_t t = (int64_t)i - 45000;

      // Convert metres to centimetres
      int64_t h = (int64_t)(altitude * 100.0);

      state.n++;

      state.sumT  += t;
      state.sumH  += h;
      state.sumTH += t * h;
      state.sumT2 += t * t;
    }

    uint32_t end = ESP.getCycleCount();

    uint32_t cycles = end - start;

    if (cycles > state.worstCycles)
      state.worstCycles = cycles;
  }


  // ===== LEAST-SQUARES CALCULATION =====

  double n    = (double)state.n;
  double sumT  = (double)state.sumT;
  double sumH  = (double)state.sumH;
  double sumTH = (double)state.sumTH;
  double sumT2 = (double)state.sumT2;

  double denominator =
      n * sumT2 - sumT * sumT;

  double slope =
      (n * sumTH - sumT * sumH) / denominator;

  double intercept =
      (sumH - slope * sumT) / n;


  // Convert cm/ms → m/s
  state.descent = -slope * 10.0;


  // t' = 0 is 45 seconds
  double landingShifted = -intercept / slope;

  // Convert back to original time
  state.landTime =
      45.0 + landingShifted / 1000.0;
}


// ===== COMMAND HANDLER =====

void handleCommand(String command) {

  command.trim();

  if (command == "DESCENT") {

    uint32_t start = ESP.getCycleCount();

    double answer = state.descent;

    uint32_t cycles =
        ESP.getCycleCount() - start;

    Serial.print("DESCENT = ");
    Serial.print(answer, 4);
    Serial.print(" m/s, cycles = ");
    Serial.println(cycles);
  }

  else if (command == "LAND") {

    uint32_t start = ESP.getCycleCount();

    double answer = state.landTime;

    uint32_t cycles =
        ESP.getCycleCount() - start;

    Serial.print("LAND = ");
    Serial.print(answer, 4);
    Serial.print(" s, cycles = ");
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

  Serial.println("Catch the Rocket - Level 1");

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