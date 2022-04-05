//PDEs and Integration
//CSCI 5611 Swinging Rope [Exercise]
//Stephen J. Guy <sjguy@umn.edu>

//Create Window
String windowTitle = "Swinging Rope";
Camera camera;
PImage img;
void setup() {
  size(1000, 1000, P3D);
  img = loadImage("grid.jpg");
  camera = new Camera();
  surface.setTitle(windowTitle);
  initScene();
}

//Simulation Parameters
float floor = 1000;
Vec3 gravity = new Vec3(0,50,0);
float radius = 2;


float obstacleRadius = 50;
Vec3 obstaclePos = new Vec3(0,95,140);


Vec3 stringTop = new Vec3(0,0,0);
float restLen = 3;
float strokeWidth = 1;
float mass = 10;
float k = 1000;
float kv = 750; 

//Initial positions and velocities of masses
static int maxNodes = 3000;
static int maxRopes = 1000;

//Multiple Ropes
Vec3 pos[][] = new Vec3[maxRopes][maxNodes];
Vec3 vel[][] = new Vec3[maxRopes][maxNodes];
Vec3 acc[][] = new Vec3[maxRopes][maxNodes];

//when 15 or more nodes, lines tend to stack
//if in debug mode --> lines are appearing pink instead of red 255
int numRopes = 50;
int numNodes = 47;







void initScene() {
  float ballSpread = restLen;
  stringTop.x -= (numRopes/2) * ballSpread;
  stringTop.x += ballSpread/2;
  for (int rope = 0; rope < numRopes; rope++) {
    for (int node = 0; node < numNodes; node++) {
      pos[rope][node] = new Vec3(0, 0, 0);
      pos[rope][node].x = stringTop.x;
      //pos[rope][node].y = stringTop.y + restLen*node;
      pos[rope][node].z = stringTop.z + restLen*node;
      vel[rope][node] = new Vec3(0, 0, 0);
    }
      stringTop.x += ballSpread;
  }
}







void update(float dt){
  
  //Reset accelerations each timestep (momenum only applies to velocity)
  for(int rope = 0; rope < numRopes; rope++){
    for (int node = 0; node < numNodes; node++){
      acc[rope][node] = new Vec3(0,0,0);
      acc[rope][node].add(gravity);
    }
  }
  
  //Compute (damped) Hooke's law for each spring Horizontal
  for(int rope = 0; rope < numRopes - 1; rope++){
    for (int node = 0; node < numNodes; node++){
      Vec3 diff = pos[rope+1][node].minus(pos[rope][node]);
      float stringF = -k*(diff.length() - restLen);
      //println(stringF,diff.length(),restLen);

      Vec3 stringDir = diff.normalized();
      float projVbot = dot(vel[rope][node], stringDir);
      float projVtop = dot(vel[rope+1][node], stringDir);
      float dampF = -kv*(projVtop - projVbot);

      Vec3 force = stringDir.times(stringF+dampF);
      acc[rope][node].add(force.times(-1.0/mass));
      acc[rope+1][node].add(force.times(1.0/mass));
    }
  }
  
  //Compute (damped) Hooke's law for each spring Vertical
  for(int rope = 0; rope < numRopes; rope++){
    for (int node = 0; node < numNodes-1; node++){
      Vec3 diff = pos[rope][node+1].minus(pos[rope][node]);
      float stringF = -k*(diff.length() - restLen);
      //println(stringF,diff.length(),restLen);
      
      Vec3 stringDir = diff.normalized();
      float projVbot = dot(vel[rope][node], stringDir);
      float projVtop = dot(vel[rope][node+1], stringDir);
      float dampF = -kv*(projVtop - projVbot);
      
      Vec3 force = stringDir.times(stringF+dampF);
      acc[rope][node].add(force.times(-1.0/mass));
      acc[rope][node+1].add(force.times(1.0/mass));
    }
  }
  /*
  // ropes diagonal down right
  for (int rope = 0; rope < numRopes-1; rope++) {
    for (int node = 0; node < numNodes-1; node++) {
      Vec3 diff = pos[rope][node].minus(pos[rope+1][node+1]);
      float stringF = -k*(diff.length() - restLen);
      //println(stringF,diff.length(),restLen);

      Vec3 stringDir = diff.normalized();
      float projVbot = dot(vel[rope][node], stringDir);
      float projVtop = dot(vel[rope+1][node+1], stringDir);
      float dampF = -kv*(projVtop - projVbot);

      Vec3 force = stringDir.times(stringF+dampF);
      acc[rope][node].add(force.times(-1.0/mass));
      acc[rope+1][node+1].add(force.times(1.0/mass));
    }
  }
  // ropes diagonal down left
  for (int rope = 1; rope < numRopes; rope++) {
    for (int node = 1; node < numNodes; node++) {
      Vec3 diff = pos[rope-1][node].minus(pos[rope][node-1]);
      float stringF = -k*(diff.length() - restLen);
      //println(stringF,diff.length(),restLen);

      Vec3 stringDir = diff.normalized();
      float projVbot = dot(vel[rope-1][node], stringDir);
      float projVtop = dot(vel[rope][node-1], stringDir);
      float dampF = -kv*(projVtop - projVbot);

      Vec3 force = stringDir.times(stringF+dampF);
      acc[rope-1][node].add(force.times(-1.0/mass));
      acc[rope][node-1].add(force.times(1.0/mass));
    }
  }
  */
  /*
  f_aero = -0.5 * den * c * |v|^2 * A * n
  
  den = density of fluid (air)
  v = velocity (relative to fluid)
  A = area
  c = drag coefficient
  f = force
  n = normal vector
  
  
  */

  //rK4 integration
  for(int rope = 0; rope < numRopes; rope++){
      for (int node = 1; node < numNodes; node++){
        Vec3 vK1 = vel[rope][node];
        Vec3 vK2 = vel[rope][node].plus(acc[rope][node].times(dt/2));
        Vec3 vK3 = vel[rope][node].plus(acc[rope][node].times(dt/2));
        Vec3 vK4 = vel[rope][node].plus(acc[rope][node].times(dt));
        
        Vec3 vK = vK1.plus((vK2.times(2)).plus((vK3.times(2)).plus(vK4)));

        pos[rope][node].add(vK.times(dt/6));
        vel[rope][node].add(acc[rope][node].times(dt));
      }
  }

/*
  //Midpoint integration
  Vec3 vHalf;
  for(int rope = 0; rope < numRopes; rope++){
    for(int node = 1; node < numNodes; node++){
      vHalf = vel[rope][node].plus(acc[rope][node].times(0.5*dt));
      pos[rope][node].add(vHalf.times(dt));
      vel[rope][node].add(acc[rope][node].times(dt));
    }
  }

  //Eulerian integration
  for(int rope = 0; rope < numRopes; rope++){
      for (int node = 1; node < numNodes; node++){
      vel[rope][node].add(acc[rope][node].times(dt));
      pos[rope][node].add(vel[rope][node].times(dt));
    }
  }
*/

  //Collision detection and response
  for (int rope = 0; rope < numRopes; rope++) {
    for (int node = 0; node < numNodes; node++) {
      if (pos[rope][node].y+radius > floor) {
        vel[rope][node].y *= -.9;
        pos[rope][node].y = floor - radius;
      }
      if (pos[rope][node].distanceTo(obstaclePos) < (obstacleRadius + radius)){
        Vec3 normal = (pos[rope][node].minus(obstaclePos)).normalized();
        pos[rope][node] = obstaclePos.plus(normal.times(obstacleRadius+radius).times(1.01));
        Vec3 velNormal = normal.times(dot(vel[rope][node],normal));
        vel[rope][node].subtract(velNormal.times(1.8));
      }
      
      //if (pos[rope][node].distanceTo(obstaclePos)+1 < (obstacleRadius+radius)) {
      //  Vec3 normal = (pos[rope][node].minus(obstaclePos)).normalized();
      //  pos[rope][node] = obstaclePos.plus(normal.times(obstacleRadius+radius).times(1.01));
      //  Vec3 velNormal = normal.times(dot(vel[rope][node], normal));
      //  vel[rope][node].subtract(velNormal.times(1.8));
      //}
    }
  }
  
}

//Draw the scene: one sphere per mass, one line connecting each pair
boolean paused = true;
boolean debug = false;
boolean move_forward = false;
boolean move_backward = false;
boolean move_left = false;
boolean move_right = false;
boolean move_up = false;
boolean move_down = false;

void draw() {
  background(235,235,235);
  noLights();
  if(debug){
    stroke(255,0,0,50);
    strokeWeight(1);
  }
  if(!debug){
    stroke(0,0,0,0);
    strokeWeight(0);
  }
  camera.Update(1.0/frameRate);  
  if (!paused){ update(0.00001*frameRate);}
  
  if(keyPressed){
    //obstacle movement
    if (move_forward) {
      obstaclePos.z -= 1;
    }
    if (move_backward) {
      obstaclePos.z += 1;
    }
    if (move_left) {
      obstaclePos.x -= 1;
    }
    if (move_right) {
      obstaclePos.x += 1;
    }
    if (move_down) {
      obstaclePos.y -= 1;
    }
    if (move_up) {
      obstaclePos.y += 1;
    }
  }
/*
  for (int rope = 0; rope < numRopes; rope++) {
    for (int node = 0; node < numNodes; node++) {
      pushMatrix();
      fill(0, node * 30, 255);
      translate(pos[rope][node].x, pos[rope][node].y, pos[rope][node].z);
      //sphere(radius);
      popMatrix();
    }
  }*/
  

  if(debug){
  //Draw Ropes
    for (int rope = 0; rope < numRopes-1; rope++) {
      for (int node = 0; node < numNodes-1; node++) {
        line(pos[rope][node].x, pos[rope][node].y, pos[rope][node].z, pos[rope+1][node].x, pos[rope+1][node].y, pos[rope+1][node].z);
        line(pos[rope][node].x, pos[rope][node].y, pos[rope][node].z, pos[rope][node+1].x, pos[rope][node+1].y, pos[rope][node+1].z);
        line(pos[rope][node].x, pos[rope][node].y, pos[rope][node].z, pos[rope+1][node+1].x, pos[rope+1][node+1].y, pos[rope+1][node+1].z);
        
        if(rope > 0){
          line(pos[rope][node].x, pos[rope][node].y, pos[rope][node].z, pos[rope-1][node+1].x, pos[rope-1][node+1].y, pos[rope-1][node+1].z);
        }
        if(rope == numRopes - 2){
          line(pos[numRopes-1][node].x, pos[numRopes-1][node].y, pos[numRopes-1][node].z, pos[numRopes-1][node+1].x, pos[numRopes-1][node+1].y, pos[numRopes-1][node+1].z);
          line(pos[numRopes-1][node].x, pos[numRopes-1][node].y, pos[numRopes-1][node].z, pos[numRopes-2][node+1].x, pos[numRopes-2][node+1].y, pos[numRopes-2][node+1].z);
        }
      }
      line(pos[rope][numNodes-1].x, pos[rope][numNodes-1].y, pos[rope][numNodes-1].z, pos[rope+1][numNodes-1].x, pos[rope+1][numNodes-1].y, pos[rope+1][numNodes-1].z);
    }
  }
  

  if(!debug){
    for (int rope = 1; rope < numRopes; rope++) {
      for (int node = 1; node < numNodes; node++) {
        beginShape();
        if((node%2 == 1 && rope%2 == 1) || (node%2 == 0 && rope%2 == 0)){fill(0, 130, 6);}
        else{fill(0,200,0);}
        vertex(pos[rope-1][node-1].x, pos[rope-1][node-1].y, pos[rope-1][node-1].z);
        vertex(pos[rope-1][node].x, pos[rope-1][node].y, pos[rope-1][node].z);
        vertex(pos[rope][node].x, pos[rope][node].y, pos[rope][node].z);
        vertex(pos[rope][node-1].x, pos[rope][node-1].y, pos[rope][node-1].z);
        endShape();
      }
    }
  }
  
  //Obstacle
  fill(0, 0, 200);
  pushMatrix();
  translate(obstaclePos.x, obstaclePos.y, obstaclePos.z);
  sphere(obstacleRadius);
  popMatrix();


  if (paused)
    surface.setTitle(windowTitle + " [PAUSED]");
  else
    surface.setTitle(windowTitle + " "+ nf(frameRate, 0, 2) + "FPS");

  //floor
  pushMatrix();
  translate(0, floor, 0);
  beginShape();
  texture(img);
  vertex(-10000, 0, -10000, 0, 0);
  vertex(10000, 0, -10000, img.width, 0);
  vertex(10000, 0, 10000, img.width, img.height);
  vertex(-10000, 0, 10000, 0, img.height);
  endShape();
  //box(100000, 1, 100000);
  popMatrix();
  
  
  
}







void keyPressed() {
    if (key == ' ') {
      paused = !paused;
    }
    if (key == 'f') {
      debug = !debug;
    }
    
    //obstacle movement
    if (key == 'i') {
      move_forward = true;
    }
    if (key == 'k') {
       move_backward = true;
    }
    if (key == 'j') {
       move_left = true;
    }
    if (key == 'l') {
       move_right = true;
    }
    if (key == 'u') {
       move_down = true;
    }
    if (key == 'o') {
       move_up = true;
    }
    
  camera.HandleKeyPressed();
}

void keyReleased()
{
    //obstacle movement
    if (key == 'i') {
      move_forward = false;
    }
    if (key == 'k') {
       move_backward = false;
    }
    if (key == 'j') {
       move_left = false;
    }
    if (key == 'l') {
       move_right = false;
    }
    if (key == 'u') {
       move_down = false;
    }
    if (key == 'o') {
       move_up = false;
    }
    
  camera.HandleKeyReleased();
}
