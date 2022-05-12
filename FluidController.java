final int N = 128;
final int iter = 4;
final int SCALE = 5;
Fluid fluid;

void settings()
{
  size(N * SCALE,N * SCALE);
}

void setup()
{
  fluid = new Fluid(0.01, 0.001 ,0.005);//water coloury
  //fluid = new Fluid(0.01, 0.002, 0.001);//wispy, kinda of like cutting flowing water, then a waterfall appears from the top
  //fluid = new Fluid(0.1, 0.0001, 0.001);//smoke trail
}

void mouseDragged()
{
  fluid.addDensity(mouseX/SCALE,mouseY/SCALE,5000); //water coloury
  //fluid.addDensity(mouseX/SCALE,mouseY/SCALE,3000);//wispy, kinda of like cutting flowing water, then a waterfall appears from the top
  //fluid.addDensity(mouseX/SCALE,mouseY/SCALE,10000);//smoke trail
  float amtX = mouseX - pmouseX;
  float amtY = mouseY - pmouseY;
  fluid.addVelocity(mouseX/SCALE,mouseY/SCALE, amtX/5, amtY/5);
}

void draw()
{
  background(0);
  //fluid.addDensity(int(0.5*width/SCALE),int(0.5*height/SCALE),10000);
  fluid.step();
  fluid.renderD();
  //fluid.fadeD();
}