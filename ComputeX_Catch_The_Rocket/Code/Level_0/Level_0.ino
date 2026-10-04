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

// Our remembered data
struct State {
  uint32_t count;
  uint32_t worstCycles;
};

State state = {0, 0};

void setup() {
  Serial.begin(115200);

  Serial.println("Catch the Rocket - Level 0");

  Serial.print("State size = ");
  Serial.print(sizeof(State));
  Serial.println(" bytes");

  for (uint32_t i = 0; i < 100000; i++) {

    uint32_t start = ESP.getCycleCount();

    float altitude = nextAltitude();

    uint32_t end = ESP.getCycleCount();

    uint32_t cycles = end - start;

    if (cycles > state.worstCycles)
      state.worstCycles = cycles;

    state.count++;

    // Prevent compiler from treating the result as unused
    if (altitude < -1000000.0f)
      Serial.println("Impossible");
  }

  Serial.print("Readings = ");
  Serial.println(state.count);

  Serial.print("Worst cycles = ");
  Serial.println(state.worstCycles);

  Serial.println("Level 0 complete.");
}

void loop() {
}