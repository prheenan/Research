int incomingByte = 0;

void setup()
{
  // Open serial connection.
  Serial.begin(9600);
}
 
void loop()
{
  if (Serial.available() > 0) 
  {
   float a = Serial.parseFloat();
   // See also: https://www.arduino.cc/en/Tutorial/StringToIntExample
   // say what you got, plus 1:
   Serial.print("Adding 1 gives: "); // ASCII printable characters
   Serial.println( (float) (a+1), 1);
  }
} 
 
