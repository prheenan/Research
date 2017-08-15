int incomingByte = 0;

void setup()
{
  // Open serial connection.
  Serial.begin(9600);
}
 
void loop()
{
  if (Serial.available() > 0) {
   // read the incoming byte:
   incomingByte = Serial.read();
   // See also: https://www.arduino.cc/en/Tutorial/StringToIntExample
   // say what you got, plus 1:
   Serial.print("Adding 1 gives: "); // ASCII printable characters
   Serial.println(incomingByte+1, DEC);
  }
} 
 
