#include "t3spi.h"
#include "hardware.h"

//Initialize T3SPI class as SPI_MASTER



//IO for boartd selection


//calibration data
float zeroPoint[NUM_B][NUM_C];
bool calibrationStatus[NUM_B][NUM_C];
bool should_next = false;
float gain[NUM_B][NUM_C];

volatile int counter = 0;
volatile int cint = 0;

int8_t channels_used[NUM_B][NUM_C]  =
{ 
  {0},
};

int8_t channel_order[NUM_B * NUM_C];
int8_t board_order[NUM_B * NUM_C];

 
const int loopsize = 1; 
float output_currents[loopsize][1] = 
{ 
{.1},
};

/******************************************************/
/*********************** OUTPUT COMPUTATION     *******/
/******************************************************/

uint16_t computeDacVal_V(float voltage, int b, int c) {
  return uint16_t((65535.0 * (voltage - zeroPoint[b][c]) / 5.0));
}

uint16_t computeDacVal_I(float current, int b, int c) {
  return uint16_t((65535.0 * (current / gain[b][c] + 2.5 - zeroPoint[b][c]) / 5.0));
}

float computeOutV(uint16_t dacVal) {
  return ((float(dacVal) * 4.096 / 4096.0));
}

float computeOutI(uint16_t dacVal) {
  return ((float(dacVal) * 4.096 / 4096.0) - 1.25) / 10 / 0.2;
}

int t;



/******************************************************/
/*********************** UTILITY *********************/
/******************************************************/
float zero_all() {
  for (int b = 0; b < NUM_B; b++) {
    selectBoard(b);
    Serial.println("---------------");
    for (int c = 0; c < NUM_C; c++) {
      LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_I(0, b, c));
      Serial.println("-------done--------");
    }
  }
}

float measure_gain(uint8_t b, uint8_t c) {
  //jump to 2.0 first so output returns nutral;

  selectBoard(b);
  delay(1);
  LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_V(2.0, 0, 0));
  delayMicroseconds(1000);
  uint16_t out_2v0 = LTC1863ReadSlow(c, 50);
  //  Serial.println(computeOutI(out_2v0),5);
  LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_V(2.5, 0, 0));
  delayMicroseconds(1000);
  uint16_t out_2v5 = LTC1863ReadSlow(c, 50);
  //  Serial.println(computeOutI(out_2v5),5);

  return (computeOutI(out_2v5) - computeOutI(out_2v0)) / (0.5);
}

bool calibrate_channel(uint8_t b, uint8_t c) {
  zeroPoint[b][c] = 0;
  delay(1);
  gain[b][c] = measure_gain(b, c);
  if (abs(gain[b][c] + 1.62) > 0.5) {
    //    gain[b][c] = 1.6;d
    Serial.println("failed (gain)");
    calibrationStatus[b][c] = false;
    return false;
  } else {
    //    return true;
  }

  //  Serial.print("gain: ");
  //  Serial.println(gain[b][c]);
  for (int i = 0; i < 10; i++) {
    float output_offset_I = computeOutI(LTC1863ReadSlow(c));
    //    Serial.print("iteration: ");
    //    Serial.println(i);
    //    Serial.println(output_offset_I,5);
    if (abs(output_offset_I) <= 0.001) {
      calibrationStatus[b][c] = true;
      return true;
    }
    zeroPoint[b][c] = zeroPoint[b][c] + (output_offset_I / gain[b][c]);
    //    Serial.print("next: ");
    //    Serial.println(zeroPoint[b][c],5);
    LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_I(0, b, c));
    delay(10);
  }
    Serial.println("failed (cal)");
  calibrationStatus[b][c] = false;
  zeroPoint[b][c] = 0;
  LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_I(0, b, c));
  return false;
}

bool calibrate_all() {
  for (int b = 0; b < NUM_B; b++) {
    for (int c = 0; c < NUM_C; c++) {
      calibrate_channel(b, c);
    }
    delay(500);
  }
}

void print_all_boards() {
  for (int b = 0; b < NUM_B; b++) {
    selectBoard(b);
    Serial.println("---------------");
    Serial.print("B: ");
    Serial.println(b);
    for (int c = 0; c < NUM_C; c++) {
      Serial.print(c);
      Serial.print(": ");
      uint16_t data = LTC1863ReadSlow(c,50);
      Serial.print(computeOutI(data), 4);
      Serial.print("\t");
      Serial.print(gain[b][c]);
      if (!calibrationStatus[b][c]) {
        Serial.println(" X");
      } else {
        Serial.println("");
      }

    }
  }
}

void update_outputs(int idx) {

  int8_t b = 0;
  selectBoard(b);
  for (int i = 0; i < NUM_C * NUM_B; i++) {

    int c = channel_order[i];
    if (c == -1 || c == -1) {
      break;
    }
    if (b != board_order[i]) {
      b = board_order[i];
      selectBoard(b);
    }
    //    Serial.println(output_currents[idx][i]);
    LTC2656Write(WRITE_AND_UPDATE, channelMap[c], computeDacVal_I(output_currents[idx][i], b, c));
  }
}

void print_all() {
  int8_t b = 0;
  selectBoard(b);
  Serial.println("-------------");
  for (int i = 0; i < NUM_C * NUM_B; i++) {

    int c = channel_order[i];
    if (c == -1) {
      break;
    }
    if (b != board_order[i]) {
      b = board_order[i];
      selectBoard(b);
    }
    uint16_t data = LTC1863ReadSlow(c);
    Serial.print(i);
    Serial.print("(");
    Serial.print(b);
    Serial.print(",");
    Serial.print(c);
    Serial.print(")\t");
    Serial.print(computeOutI(data), 4);
    Serial.print("\t");
    Serial.print(gain[b][c]);
    if (!calibrationStatus[b][c]) {
      Serial.println(" X");
    } else {
      Serial.println("");
    }
  }
}

/******************************************************/
/*********************** IRUPT ************************/
/******************************************************/
void setDACVal() {
  update_outputs(counter);
  Serial.println(counter);
  counter+=1;

  if (counter%loopsize== 0 ) {
    counter = 0;
  }
//  if (counter <= 0) {
//    counter = 39;
//  }
}

//}


/******************************************************/
/*********************** SETUP ************************/
/******************************************************/

void setup() {
  read_in_flight = false;
  Serial.begin(115200);

  //SETUP board and function select
  initIO();
  selectNone();
  spiInit();

  //Initialize calibration data
  //for (int b = 0; b < NUM_B; b++) {
    selectBoard(0);
    for (int c = 0; c < NUM_C; c++) {
      zeroPoint[2][c] = 0;
      gain[2][c] = -1.6;
      //      LTC2656Write(WRITE_AND_UPDATE,channelMap[c],32768);
    }
  //}


  for (int j = 0; j < NUM_B * NUM_C; j++) {
    channel_order[j] = -1;
    board_order[j] = -1;
  }
  int i = 0;
  for (int b = 0; b < NUM_B; b++) {
    for (int c = 0; c < NUM_C; c++) {
      if (channels_used[b][c] != -1) {
        channel_order[i] = channels_used[b][c];
        board_order[i] = b;
        i = i + 1;
      } else {
        break;
      }
    }
  }

  delay(500);
//    attachInterrupt(interruptPin, setDACVal, FALLING);
//    attachInterrupt(4, setDACVal, RISING);
  Serial.println("I'm up");
  selectBoard(0);

  delay(100);
  NVIC_ENABLE_IRQ(IRQ_SPI1);
  t = micros();
  Serial.println("about to write");
  //  zero_all();
  Serial.println("finihsed write");
  delay(500);
  SPI_SLAVE->packetCT = 0;
  SPI_SLAVE->dataPointer = 0;

  //Serial.println(measure_gain(2,0),1);

  Serial.println("lets do it");
  selectBoard(0);
  Serial.println("lets do it2");
//  zero_all();
//  Serial.println("lets do it3");
//  selectBoard(0);
  Serial.println("lets do it4");
  delay(100);
  Serial.println(NUM_B);
  Serial.println(NUM_C);
}


/******************************************************/
/*********************** MAIN LOOP ********************/
/******************************************************/

void loop() {
  //Serial.println("looping");

  if (should_next) {
    //    zero_all();
    //    calibrate_all();
    //    update_outputs(counter);
    should_next = false;
    print_all();
  }

  char incomingByte;

  if (Serial.available() > 0) {
    Serial.println("got a byte");

    incomingByte = Serial.read();
    Serial.print(incomingByte);
    switch (incomingByte) {
      case 'Z':
        zero_all();
        break;
      case 'C':
        calibrate_channel(2, 0);
        //for (int i = 0; i < 8; i++) {
        //  calibrate_channel(0, i);
        //}
        //for (int i = 0; i < 8; i++) {
        //  calibrate_channel(1, i);
        //}
        //for (int i = 0; i < 8; i++) {
        //  calibrate_channel(2, i);
        //}
        //for (int i = 0; i < 8; i++) {
        //  calibrate_channel(3, i);
        //}
        break;
      case 'D':
        for (int i = 0; i < 8; i++) {
          calibrate_channel(0, i);
        }
        break;
      case 'I':
        print_all();
        break;
      case 'A':
        print_all_boards();
        break;
      case 'S':
        selectBoard(0);
        LTC2656Write(WRITE_AND_UPDATE, channelMap[0], computeDacVal_I(0.5, 0, 0));
        break;
      case 'M':
        counter = 0;
        cint = 0;
        should_next = false;
        Serial.println(counter);
        update_outputs(counter);
        counter+=1;

    }
  }


}
