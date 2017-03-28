/*
 *	Includes
 */
#include <algorithm>
#include <assimp/cimport.h>		// C importer
#include <assimp/scene.h>		// Collects data
#include <assimp/postprocess.h> // Various extra operations
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <iomanip>				// setprecision
#include <iostream>
#include <math.h>
#include <mmsystem.h>
#include <queue>
#include <sstream>
#include <stdio.h>
#include <vector>				// STL dynamic memory
#include <windows.h>

#include "Antons_maths_funcs.h" // Anton's maths functions
#include "Camera.h"
#include "Distance.h"
#include "Mesh.h"
#include "Model.h"
#include "RigidBody.h"
#include "Shader_Functions.h"
#include "text.h"
#include "time.h"

using namespace std;

/*
 *	Globally defined variables and constants
 */
#define NUM_MESHES   8
#define NUM_SHADERS	 5
#define NUM_TEXTURES 7

bool firstMouse = true;
bool keys[1024];
bool pause = false;
Camera camera(vec3(0.0f, 0.0f, 4.0f));
char name1[20], name2[20];
const GLuint numRigidBodies = 2;
enum Meshes { PLANE_MESH, D6_MESH, SPHERE_MESH, D4_MESH, D8_MESH, D10_MESH, D12_MESH, D20_MESH};
enum Shaders { SKYBOX, BASIC_COLOUR_SHADER, BASIC_TEXTURE_SHADER, LIGHT_SHADER, LIGHT_TEXTURE_SHADER };
enum Textures { PLANE_TEXTURE, D4_TEXTURE, D6_TEXTURE, D8_TEXTURE, D10_TEXTURE, D12_TEXTURE, D20_TEXTURE};
GLfloat cameraSpeed = 0.005f;
GLfloat currentDistance = -1.0f;
GLuint lastX = 400, lastY = 300;
GLuint mode = AABB;
GLuint shaderProgramID[NUM_SHADERS];
int screenWidth = 1000;
int screenHeight = 800;
int stringIDs[5];
Mesh asteroid, boundingBox, sphereMesh, pyramidMesh, d4, d6, d8, d10, d12, d20;
Mesh vertexMesh, edgeMesh, faceMesh;
vec4 backgroundColours[3] = { vec4(0.0f, 0.0f, 0.0f, 1.0f), vec4(0.1f, 0.0f, 0.0f, 1.0f), vec4(0.2f, 0.0f, 0.0f, 1.0f) };
vec4 featureColours[2] = { vec4(1.0f, 1.0f, 0.0f, 1.0f), vec4(0.0f, 1.0f, 1.0f, 1.0f) };
vector<RigidBody> rigidbodies;
void *feature1, *feature2;
void *features[2];

// | Resource Locations
const char * meshFiles[NUM_MESHES] = { "../Meshes/plane.obj", "../Meshes/cube.dae", "../Meshes/particle.dae", "../Meshes/d4.obj", "../Meshes/d8.obj", "../Meshes/d10.obj", "../Meshes/d12.obj", "../Meshes/d20.obj" };
const char * skyboxTextureFiles[6] = { "../Textures/DSposx.png", "../Textures/DSnegx.png", "../Textures/DSposy.png", "../Textures/DSnegy.png", "../Textures/DSposz.png", "../Textures/DSnegz.png"};
const char * textureFiles[NUM_TEXTURES] = { "../Textures/plane.jpg", "../Textures/d4.png", "../Textures/d6.png", "../Textures/d8.png", "../Textures/d10.png", "../Textures/d12.png", "../Textures/d20.png" };

const char * vertexShaderNames[NUM_SHADERS] = { "../Shaders/SkyboxVertexShader.txt", "../Shaders/ParticleVertexShader.txt", "../Shaders/BasicTextureVertexShader.txt", "../Shaders/LightVertexShader.txt", "../Shaders/LightTextureVertexShader.txt" };
const char * fragmentShaderNames[NUM_SHADERS] = { "../Shaders/SkyboxFragmentShader.txt", "../Shaders/ParticleFragmentShader.txt", "../Shaders/BasicTextureFragmentShader.txt", "../Shaders/LightFragmentShader.txt", "../Shaders/LightTextureFragmentShader.txt" };

void initialiseRigidBodies(bool indices);

string frf(const float &f)
{
	ostringstream ss;
	ss << setfill(' ') << std::setw(6) << fixed << setprecision(3) << f;
	string s(ss.str());
	return s;
}

void draw_text()
{
	ostringstream typeOSS, numOSS, distanceOSS, feature1OSS, feature2OSS;
	string typeString, numString, distanceString, feature1String, feature2String;
	if (mode == BOUNDING_SPHERES)
		typeOSS << "Broad Phase Collision: Bounding Spheres";
	else if (mode == AABB)
		typeOSS << "Broad Phase Collision: AABB";
	else
		typeOSS << "Broad Phase Collision: None";

	numOSS << "Number of rigid bodies: " << numRigidBodies;
	if (currentDistance < 0.0f)
		distanceOSS << "Distance between meshes: NA";
	else
		distanceOSS << "Distance between meshes: " << currentDistance;

	feature1OSS << "Feature 1: " << name1;
	feature2OSS << "Feature 2: " << name2;
	
	typeString = typeOSS.str();
	numString = numOSS.str();
	distanceString = distanceOSS.str();
	feature1String = feature1OSS.str();
	feature2String = feature2OSS.str();
	
	update_text(stringIDs[0], typeString.c_str());
	update_text(stringIDs[1], numString.c_str());
	update_text(stringIDs[2], distanceString.c_str());
	update_text(stringIDs[3], feature1String.c_str());
	update_text(stringIDs[4], feature2String.c_str());

	draw_texts();
}

void init_text()
{
	stringIDs[0] = add_text("Broad Phase Collision: ", -0.95f, 0.95f, 25.0f, 1.0f, 1.0f, 1.0f, 1.0f);
	stringIDs[1] = add_text("Number of rigid bodies: ", -0.95f, 0.9f, 25.0f, 1.0f, 1.0f, 1.0f, 1.0f);
	stringIDs[2] = add_text("Distance between meshes: NA", -0.95f, 0.85f, 25.0f, 1.0f, 1.0f, 1.0f, 1.0f);
	stringIDs[3] = add_text("Feature 1: NA", -0.95f, 0.8f, 25.0f, 1.0f, 1.0f, 1.0f, 1.0f);
	stringIDs[4] = add_text("Feature 2: NA", -0.95f, 0.75f, 25.0f, 1.0f, 1.0f, 1.0f, 1.0f);
}

void drawClosestFeatures(mat4 view, mat4 projection)
{
	for (int i = 0; i < 2; i++)
	{
		int featureType = featureTag(features[i]);

		if (featureType == V)
		{
			vec4 p1 = vec4(((vertex *)features[i])->coords, 0.0f);
			p1 = (rigidbodies[i].rotation * p1) + rigidbodies[i].position;
			GLfloat vertex_array[] = {
				p1.v[0], p1.v[1], p1.v[2]
			};
			vertexMesh = Mesh(&shaderProgramID[BASIC_COLOUR_SHADER]);
			vertexMesh.generateObjectBufferMesh(vertex_array, 1);

			mat4 vertex_model = identity_mat4();
			glPointSize(20.0f);
			vertexMesh.drawPoint(view, projection, vertex_model, featureColours[i]);
			glPointSize(1.0f);
		}
		else if (featureType == E)
		{
			vec4 p1 = vec4(((edge *)features[i])->v1->coords, 0.0f);
			p1 = (rigidbodies[i].rotation * p1) + rigidbodies[i].position;
			vec4 p2 = vec4(((edge *)features[i])->v2->coords, 0.0f);
			p2 = (rigidbodies[i].rotation * p2) + rigidbodies[i].position;

			GLfloat edge_array[] = {
				p1.v[0], p1.v[1], p1.v[2],
				p2.v[0], p2.v[1], p2.v[2]
			};
			edgeMesh = Mesh(&shaderProgramID[BASIC_COLOUR_SHADER]);
			edgeMesh.generateObjectBufferMesh(edge_array, 2);

			mat4 edge_model = identity_mat4();
			glLineWidth(10.0f);
			edgeMesh.drawLine(view, projection, edge_model, featureColours[i]);
			glLineWidth(1.0f);
		}
		else if (featureType == F)
		{
			vec4 p1 = vec4(((face *)features[i])->verts->f.v->coords, 0.0f);
			p1 = (rigidbodies[i].rotation * p1) + rigidbodies[i].position;
			vec4 p2 = vec4(((face *)features[i])->verts->next->f.v->coords, 0.0f);
			p2 = (rigidbodies[i].rotation * p2) + rigidbodies[i].position;
			vec4 p3 = vec4(((face *)features[i])->verts->next->next->f.v->coords, 0.0f);
			p3 = (rigidbodies[i].rotation * p3) + rigidbodies[i].position;

			GLfloat face_array[] = {
				p1.v[0], p1.v[1], p1.v[2],
				p2.v[0], p2.v[1], p2.v[2],
				p3.v[0], p3.v[1], p3.v[2]
			};
			faceMesh = Mesh(&shaderProgramID[BASIC_COLOUR_SHADER]);
			faceMesh.generateObjectBufferMesh(face_array, 3);

			mat4 face_model = identity_mat4();
			//face_model = scale(face_model, vec3(1.005f, 1.005f, 1.005f));
			faceMesh.drawMesh(view, projection, face_model, featureColours[i]);
		}
	}
}

void display() 
{
	// Tell GL to only draw onto a pixel if the shape is closer to the viewer
	glEnable(GL_DEPTH_TEST);	// Enable depth-testing
	glDepthFunc(GL_LESS);		// Depth-testing interprets a smaller value as "closer"
	if (currentDistance < 0.005f && currentDistance > 0.0f)
		glClearColor(backgroundColours[2].v[0], 0.0f, 0.0f, 1.0f);
	else if(currentDistance < 0.01f && currentDistance > 0.0f)
		glClearColor(backgroundColours[1].v[0], 0.0f, 0.0f, 1.0f);
	else
		glClearColor(backgroundColours[0].v[0], 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	mat4 view = camera.GetViewMatrix();
	mat4 projection = perspective(camera.Zoom, (float)screenWidth / (float)screenHeight, 0.1f, 100.0f);
	mat4 model = identity_mat4();
	vec4 view_position = vec4(camera.Position.v[0], camera.Position.v[1], camera.Position.v[2], 0.0f);

	if (rigidbodies[0].collisionAABB)
	{
		drawClosestFeatures(view, projection);
	}

	for (GLuint i = 0; i < numRigidBodies; i++)
	{
		rigidbodies[i].drawMesh(view, projection, view_position);
		if (mode == BOUNDING_SPHERES)
			rigidbodies[i].drawBoundingSphere(view, projection);
		else if (mode == AABB)
			rigidbodies[i].drawAABB(view, projection, &shaderProgramID[BASIC_COLOUR_SHADER]);
	}
	
	boundingBox.drawLine(view, projection, model, vec4(1.0f, 1.0f, 0.0f, 1.0f));


	draw_text();
	
	glutSwapBuffers();
}

void processInput()
{
	if (keys[GLUT_KEY_UP])
		camera.ProcessKeyboard(FORWARD, cameraSpeed);
	if(keys[GLUT_KEY_DOWN])
		camera.ProcessKeyboard(BACKWARD, cameraSpeed);
	if (keys[GLUT_KEY_LEFT])
		camera.ProcessKeyboard(LEFT, cameraSpeed);
	if (keys[GLUT_KEY_RIGHT])
		camera.ProcessKeyboard(RIGHT, cameraSpeed);

	if (keys['1'])
		mode = BOUNDING_SPHERES;
	if (keys['2'])
		mode = AABB;

	if (keys['5'])
		useForce = -1;
	if (keys['6'])
		useForce = FORCE1;
	if (keys['7'])
		useForce = FORCE2;
	if (keys['8'])
		useForce = FORCE3;
	if (keys['9'])
		useForce = FORCE4;

	if (keys['p'])
		pause = true;
	if (keys['o'])
		pause = false;

	if (keys['u'])
		rigidbodies[0].position.v[1] += 0.001f;
	if (keys['j'])
		rigidbodies[0].position.v[1] -= 0.001f;
	if (keys['y'])
		rigidbodies[0].position.v[2] += 0.001f;
	if (keys['i'])
		rigidbodies[0].position.v[2] -= 0.001f;
	if (keys['k'])
		rigidbodies[0].position.v[0] += 0.001f;
	if (keys['h'])
		rigidbodies[0].position.v[0] -= 0.001f;

	if (keys[(char)27])
		exit(0);
}

void updateScene()
{
	processInput();

	if (!pause)
	{
		updateRigidBodies(mode, numRigidBodies, rigidbodies);

		if (mode == BOUNDING_SPHERES)
			checkBoundingSphereCollisions(numRigidBodies, rigidbodies);
		else if (mode == AABB)
		{
			checkAABBCollisions(numRigidBodies, rigidbodies);
			if (rigidbodies[0].collisionAABB)
			{
				rigidbodies[0].poly.pose = rigidbodies[0].transformationMatrix;
				rigidbodies[1].poly.pose = rigidbodies[1].transformationMatrix;
				if (features[0] == NULL || features[1] == NULL)
					currentDistance = closestFeaturesInit(&rigidbodies[0].poly, &features[0], &rigidbodies[1].poly, &features[1]);
				else
					currentDistance = closestFeatures(&rigidbodies[0].poly, &features[0], &rigidbodies[1].poly, &features[1]);
				featureName(features[0], name1);
				featureName(features[1], name2);
			}
			else
			{
				features[0] = NULL;
				features[1] = NULL;
				strcpy_s(name1, sizeof(name1), "NA");
				strcpy_s(name2, sizeof(name2), "NA");
				currentDistance = -1.0f;
			}
		}
	}

	// Draw the next frame
	glutPostRedisplay();
}

void initialiseRigidBodies(bool indices)
{
	for (GLuint i = 0; i < numRigidBodies; i++)
	{
		RigidBody &rigidBody = rigidbodies[i];
		GLfloat randomX1 = ((rand() % 8) - 4) / 5.0f;
		GLfloat randomY1 = ((rand() % 8)) / 10.0f;
		GLfloat randomZ1 = ((rand() % 8) - 4) / 5.0f;

		rigidBody.position = vec4(randomX1, randomY1, randomZ1, 0.0f);

		GLfloat rand0 = (rand() % 100) / 100.0f;
		GLfloat rand1 = (rand() % 100) / 100.0f;
		GLfloat rand2 = (rand() % 100) / 100.0f;
		GLfloat rand3 = (rand() % 100) / 100.0f;

		rigidBody.orientation.q[0] = rand0;
		rigidBody.orientation.q[1] = rand1;
		rigidBody.orientation.q[2] = rand2;
		rigidBody.orientation.q[3] = rand3;
		normalise(rigidBody.orientation);

		if (indices)
		{
			rigidBody.xMinI = 2 * i;
			rigidBody.xMaxI = (2 * i) + 1;
			rigidBody.yMinI = 2 * i;
			rigidBody.yMaxI = (2 * i) + 1;
			rigidBody.zMinI = 2 * i;
			rigidBody.zMaxI = (2 * i) + 1;
		}
	}
}

void init()
{
	if (!init_text_rendering("../Textures/freemono.png", "../Textures/freemono.meta", screenWidth, screenHeight))
	{
		fprintf(stderr, "ERROR init text rendering\n");
		exit(1);
	}
	init_text();

	// Compile the shaders
	for (int i = 0; i < NUM_SHADERS; i++)
	{
		shaderProgramID[i] = CompileShaders(vertexShaderNames[i], fragmentShaderNames[i]);
	}

	GLfloat bounding_box_vertices[] = {
		1.0f, 1.0f, 1.0f,
		-1.0f, 1.0f, 1.0f,
		-1.0f, -1.0f, 1.0f,
		1.0f, -1.0f, 1.0f,

		1.0f, 1.0f, 1.0f,
		1.0f, 1.0f, -1.0f,
		1.0f, -1.0f, -1.0f,
		1.0f, -1.0f, 1.0f,

		-1.0f, -1.0f, 1.0f,
		-1.0f, -1.0f, -1.0f,
		1.0f, -1.0f, -1.0f,
		1.0f, 1.0f, -1.0f,

		-1.0f, 1.0f, -1.0f,
		-1.0f, -1.0f, -1.0f,
		-1.0f, 1.0f, -1.0f,
		-1.0f, 1.0f, 1.0f
	};

	boundingBox = Mesh(&shaderProgramID[BASIC_COLOUR_SHADER]);
	boundingBox.generateObjectBufferMesh(bounding_box_vertices, 16);

	d4 = Mesh(&shaderProgramID[BASIC_TEXTURE_SHADER]);
	d4.generateObjectBufferMesh(meshFiles[D4_MESH]);
	d4.loadTexture(textureFiles[D4_TEXTURE]);

	/*d6 = Mesh(&shaderProgramID[BASIC_TEXTURE_SHADER]);
	d6.generateObjectBufferMesh(meshFiles[D6_MESH]);
	d6.loadTexture(textureFiles[D6_TEXTURE]);*/

	d8 = Mesh(&shaderProgramID[BASIC_TEXTURE_SHADER]);
	d8.generateObjectBufferMesh(meshFiles[D8_MESH]);
	d8.loadTexture(textureFiles[D8_TEXTURE]);

	/*d10 = Mesh(&shaderProgramID[BASIC_TEXTURE_SHADER]);
	d10.generateObjectBufferMesh(meshFiles[D10_MESH]);
	d10.loadTexture(textureFiles[D10_TEXTURE]);

	d12 = Mesh(&shaderProgramID[BASIC_TEXTURE_SHADER]);
	d12.generateObjectBufferMesh(meshFiles[D12_MESH]);
	d12.loadTexture(textureFiles[D12_TEXTURE]);

	d20 = Mesh(&shaderProgramID[BASIC_TEXTURE_SHADER]);
	d20.generateObjectBufferMesh(meshFiles[D20_MESH]);
	d20.loadTexture(textureFiles[D20_TEXTURE]);*/

	sphereMesh = Mesh(&shaderProgramID[BASIC_COLOUR_SHADER]);
	sphereMesh.generateObjectBufferMesh(meshFiles[SPHERE_MESH]);

	//RigidBody rigidBody = RigidBody(asteroid.vertex_count, asteroid.vertex_positions);
	RigidBody d4rb = RigidBody(d4, 0.4f);
	d4rb.addBoundingSphere(sphereMesh, green);
	rigidbodies.push_back(d4rb);

	/*RigidBody d6rb = RigidBody(d6, 0.2f);
	d6rb.addBoundingSphere(sphereMesh, green);
	rigidbodies.push_back(d6rb);*/

	RigidBody d8rb = RigidBody(d8, 0.5f);
	d8rb.addBoundingSphere(sphereMesh, green);
	rigidbodies.push_back(d8rb);

	/*RigidBody d10rb = RigidBody(d10, 0.5f);
	d10rb.addBoundingSphere(sphereMesh, green);
	rigidbodies.push_back(d10rb);

	RigidBody d20rb = RigidBody(d20, 0.4f);
	d20rb.addBoundingSphere(sphereMesh, green);
	rigidbodies.push_back(d20rb);*/

	initialiseRigidBodies(true);
}

/*
 *	User Input Functions
 */
#pragma region USER_INPUT_FUNCTIONS
void pressNormalKeys(unsigned char key, int x, int y)
{
	keys[key] = true;
	if (keys['k'])
		useGravity = !useGravity;

	if (keys['0'])
		initialiseRigidBodies(false);
}

void releaseNormalKeys(unsigned char key, int x, int y)
{
	keys[key] = false;
}

void pressSpecialKeys(int key, int x, int y)
{
	keys[key] = true;
}

void releaseSpecialKeys(int key, int x, int y)
{
	keys[key] = false;
}

void mouseClick(int button, int state, int x, int y)
{}

void processMouse(int x, int y)
{
	if (firstMouse)
	{
		lastX = x;
		lastY = y;
		firstMouse = false;
	}

	int xoffset = x - lastX;
	int yoffset = lastY - y;

	lastX = x;
	lastY = y;

	camera.ProcessMouseMovement((GLfloat)xoffset, (GLfloat)yoffset);
}

void mouseWheel(int button, int dir, int x, int y)
{}
#pragma endregion

/*
 *	Main
 */
int main(int argc, char** argv) 
{
	srand((unsigned int)time(NULL));

	// Set up the window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(screenWidth, screenHeight);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - screenWidth) / 2, (glutGet(GLUT_SCREEN_HEIGHT) - screenHeight) / 4);
	glutCreateWindow("I-Collide Implementation");

	// Glut display and update functions
	glutDisplayFunc(display);
	glutIdleFunc(updateScene);

	// User input functions
	glutKeyboardFunc(pressNormalKeys);
	glutKeyboardUpFunc(releaseNormalKeys);
	glutSpecialFunc(pressSpecialKeys);
	glutSpecialUpFunc(releaseSpecialKeys);
	glutMouseFunc(mouseClick);
	glutPassiveMotionFunc(processMouse);
	glutMouseWheelFunc(mouseWheel);


	glewExperimental = GL_TRUE; //for non-lab machines, this line gives better modern GL support
	
	// A call to glewInit() must be done after glut is initialized!
	GLenum res = glewInit();
	// Check for any errors
	if (res != GLEW_OK) {
		fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
		return 1;
	}

	// Set up meshes and shaders
	init();
	// Begin infinite event loop
	glutMainLoop();
	return 0;
}