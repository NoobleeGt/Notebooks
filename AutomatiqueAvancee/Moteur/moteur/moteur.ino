// pin definitions
const int mv_pin = 3;
const int sp_pin = 5;
const int pv_pin = A0;

// control related variables
double pv = 0;
double sp = 1;
double e = 0;
double i = 0;
double i_prec = 0;
double mv = 0;
double mv_i = 0;

// time variable
unsigned long t_prec = 0;

void setup() {
  // put your setup code here, to run once:
  pinMode(mv_pin, OUTPUT);
  pinMode(sp_pin, OUTPUT);
    
  Serial.begin(9600);
  while (!Serial){};
}

void loop() {
  // put your main code here, to run repeatedly:
  // read a new set point on the Serial
  if (Serial.available() > 0){
    sp = Serial.readString().toInt() * 5. / 1023;
  }
  
  unsigned long t = millis();
  
  // check if 10ms has elapsed since last code execution
  if (t - t_prec >= 10){
    // read pv from analog input and scale it to 0-5 V
    pv = analogRead(pv_pin) * 5. / 1023;
    
    // PI algorithm starts here
    e = sp - pv;
    
    // integration error handling
    i = i_prec + e;
    mv_i = 0.01 / 0.3 * i;
  
    mv = 3.33 * (e + mv_i);
    
    // ARW (anti reset windup) handling
    if ((mv > 5.) || (mv < 0.)){
      i = i_prec;
      mv_i = 0.01 / 0.3 * i;
  
      mv = 3.33 * (e + mv_i);
    }
    
    // limits the manipulated value to 0-5 V
    mv = constrain(mv, 0., 5.);
    
    // write mv and pv so that we can apply it and monitor it
    analogWrite(mv_pin, round(mv * 255 / 5.));
    analogWrite(sp_pin, round(sp * 255 / 5.));
    
    // update variables
    i_prec = i;
    t_prec = t;
  }
}
