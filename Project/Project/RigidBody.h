#pragma once

#include <GL/glew.h>
#include <vector>

#include "Antons_maths_funcs.h"
#include "Distance.h"
#include "Mesh.h"

// Colours
vec4 red = vec4(1.0f, 0.0f, 0.0f, 1.0f);
vec4 green = vec4(0.0f, 1.0f, 0.0f, 1.0f);

// Plane axes
vec4 yAxis = vec4(0.0f, 1.0f, 0.0f, 0.0f);
vec4 xAxis = vec4(1.0f, 0.0f, 0.0f, 0.0f);
vec4 zAxis = vec4(0.0f, 0.0f, 1.0f, 0.0f);

// Parameters
bool useGravity = true;
enum Forces { FORCE1, FORCE2, FORCE3, FORCE4 };
enum Mode { BOUNDING_SPHERES, AABB };
GLfloat dampingFactor = 0.5f;
GLfloat deltaTime = 1.0f / 300.0f;
GLfloat friction = 0.9f;
GLfloat restitution = 0.6f; //0.996f is a close approximation of perfect elasticity
GLuint useForce = -1;
vec4 gravity = vec4(0.0f, -9.81f, 0.0f, 0.0f);

#pragma region RIGIDBODY_CLASS
class RigidBody {
public:
	// Constants
	GLfloat mass;
	mat4 Ibody;
	mat4 IbodyInv;
	vec4 bodyCOM;
	vec4 worldCOM;

	// State Variables
	vec4 position;			// x(t)
	versor orientation;		// q(t)
	vec4 linearMomentum;	// P(t)
	vec4 angularMomentum;	// L(t)

	// Derived Quantities
	mat4 Iinv;				// I-1(t) = R(t) * IbodyInv * R(t)T
	mat4 rotation;			// Rotation Matrix R(t)
	vec4 velocity;			// v(t)  = P(t) / mass
	vec4 angularVelocity;	// w(t)  = I-1(t) * L(t)

	// Computed Quantities
	vec4 torque;			// T(t)
	vec4 force;				// F(t)

	// Mesh information
	GLuint numTriangles;
	GLuint numPoints;
	vector<vec4> bodyVertices;
	vector<vec4> initialWorldVertices;
	vector<vec4> worldVertices;

	Mesh rigidBodyMesh;

	// Bounding Sphere Variables
	Mesh boundingSphere;
	GLfloat boundingSphereRadius;
	vec4 boundingSphereColour;
	GLuint collidingWith;
	vec4 bodyCentroid;
	vec4 worldCentroid;

	// AABB Variables
	GLfloat xMin, xMax, yMin, yMax, zMin, zMax;
	GLuint xMinI, xMaxI, yMinI, yMaxI, zMinI, zMaxI; // Array indices
	vec4 boundingBoxColour;
	bool collisionAABB;

	GLfloat scaleFactor;

	vec4 meshColour;

	RigidBody();
	RigidBody(Mesh rigidBodyMesh, GLfloat scaleFactor);
	RigidBody(int vertex_count, vector<float> vertex_positions);
	~RigidBody();
	GLfloat calculateBoundingSphereRadius();
	vec4 getCentroid();
	void addBoundingSphere(Mesh boundingSphere, vec4 colour);
	void computeMassInertia(bool bodyCoords);
	void drawMesh(mat4 view, mat4 projection, vec4 viewPosition);
	void drawBoundingSphere(mat4 view, mat4 projection);
	void drawAABB(mat4 view, mat4 projection, GLuint* shaderID);
};

RigidBody::RigidBody()
{
	this->mass = 1.0f;

	this->orientation.q[0] = 0.0f;
	this->orientation.q[1] = 0.0f;
	this->orientation.q[2] = 1.0f;
	this->orientation.q[3] = 0.0f;

	this->position = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->rotation = identity_mat4();
	this->linearMomentum = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->angularMomentum = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	this->torque = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->force = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	this->velocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->angularVelocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->Iinv = identity_mat4();

	this->xMin = 0.0f;
	this->xMax = 0.0f;
	this->yMin = 0.0f;
	this->yMax = 0.0f;
	this->zMin = 0.0f;
	this->zMax = 0.0f;

	this->collisionAABB = false;

	this->meshColour = vec4(0.0f, 0.0f, 0.0f, 0.0f);
}

RigidBody::RigidBody(Mesh rigidBodyMesh, GLfloat scaleFactor = 1.0f)
{
	this->rigidBodyMesh = rigidBodyMesh;
	this->scaleFactor = scaleFactor;
	int vertex_count = rigidBodyMesh.vertex_count;
	vector<float> vertex_positions = rigidBodyMesh.vertex_positions;

	// Mesh Information
	this->bodyVertices.clear();
	this->initialWorldVertices.clear();
	this->worldVertices.clear();

	for (int i = 0; i < vertex_count; i++)
	{
		this->bodyVertices.push_back(vec4(vertex_positions[i * 3], vertex_positions[1 + i * 3], vertex_positions[2 + i * 3], 0.0f) * scaleFactor);
		this->initialWorldVertices.push_back(vec4(vertex_positions[i * 3], vertex_positions[1 + i * 3], vertex_positions[2 + i * 3], 0.0f) * scaleFactor);
	}

	sort(this->initialWorldVertices.begin(), this->initialWorldVertices.end());
	this->initialWorldVertices.erase(unique(this->initialWorldVertices.begin(), this->initialWorldVertices.end()), this->initialWorldVertices.end());
	this->numPoints = this->initialWorldVertices.size();
	this->worldVertices = this->initialWorldVertices;

	this->numTriangles = vertex_count / 3;

	// Constants
	computeMassInertia(false);
	this->IbodyInv = inverse(this->Ibody);

	// State Variables
	this->worldCOM = this->bodyCOM;
	this->position = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->orientation.q[0] = 0.0f;
	this->orientation.q[1] = 0.0f;
	this->orientation.q[2] = 1.0f;
	this->orientation.q[3] = 0.0f;
	this->linearMomentum = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->angularMomentum = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	// Derived Quantities
	this->rotation = quat_to_mat4(this->orientation);
	this->Iinv = this->rotation * this->IbodyInv * transpose(this->rotation);
	this->velocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->angularVelocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	// Computed Quantities
	this->torque = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->force = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	this->xMin = 0.0f;
	this->xMax = 0.0f;
	this->yMin = 0.0f;
	this->yMax = 0.0f;
	this->zMin = 0.0f;
	this->zMax = 0.0f;

	this->collisionAABB = false;

	this->meshColour = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	this->bodyCentroid = getCentroid();
}

RigidBody::RigidBody(int vertex_count, vector<float> vertex_positions)
{
	// Mesh Information
	this->bodyVertices.clear();
	this->worldVertices.clear();

	for (int i = 0; i < vertex_count; i++)
	{
		this->bodyVertices.push_back(vec4(vertex_positions[i * 3], vertex_positions[1 + i * 3], vertex_positions[2 + i * 3], 0.0f));
		this->worldVertices.push_back(vec4(vertex_positions[i * 3], vertex_positions[1 + i * 3], vertex_positions[2 + i * 3], 0.0f));
	}

	sort(this->worldVertices.begin(), this->worldVertices.end());
	this->worldVertices.erase(unique(this->worldVertices.begin(), this->worldVertices.end()), this->worldVertices.end());
	this->numPoints = this->worldVertices.size();

	this->numTriangles = vertex_count / 3;

	// Constants
	computeMassInertia(false);
	this->IbodyInv = inverse(this->Ibody);

	// State Variables
	this->worldCOM = this->bodyCOM;
	this->position = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->orientation.q[0] = 0.0f;
	this->orientation.q[1] = 0.0f;
	this->orientation.q[2] = 1.0f;
	this->orientation.q[3] = 0.0f;
	this->linearMomentum = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->angularMomentum = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	// Derived Quantities
	this->rotation = quat_to_mat4(this->orientation);
	this->Iinv = this->rotation * this->IbodyInv * transpose(this->rotation);
	this->velocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->angularVelocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	// Computed Quantities
	this->torque = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	this->force = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	this->xMin = 0.0f;
	this->xMax = 0.0f;
	this->yMin = 0.0f;
	this->yMax = 0.0f;
	this->zMin = 0.0f;
	this->zMax = 0.0f;

	this->collisionAABB = false;

	this->meshColour = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	this->bodyCentroid = getCentroid();
}

RigidBody::~RigidBody()
{}

void RigidBody::addBoundingSphere(Mesh boundingSphere, vec4 colour)
{
	this->boundingSphere = boundingSphere;
	this->boundingSphereColour = colour;
	this->boundingSphereRadius = calculateBoundingSphereRadius();
}

GLfloat RigidBody::calculateBoundingSphereRadius()
{
	GLfloat max = 0.0f;
	for (GLuint i = 0; i < this->numPoints; i++)
	{
		if (getDistance(this->worldVertices[i], this->bodyCentroid) > max)
			max = getDistance(this->worldVertices[i], this->bodyCentroid);
	}
	return max;
}

vec4 RigidBody::getCentroid()
{
	vec4 sum = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	for (GLuint i = 0; i < this->numPoints; i++)
	{
		sum += this->worldVertices[i];
	}
	return sum/(float)this->numPoints;
}

// Adapted from:
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

void RigidBody::computeMassInertia(bool bodyCoords)
{
	float oneDiv6 = (1.0f / 6.0f);
	float oneDiv24 = (1.0f / 24.0f);
	float oneDiv60 = (1.0f / 60.0f);
	float oneDiv120 = (1.0f / 120.0f);

	// order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx
	float integral[10] = { 0.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };

	int index = 0;
	for (GLuint i = 0; i < numTriangles; ++i)
	{
		// Get vertices of triangle i.
		vec3 v0 = bodyVertices[index++];
		vec3 v1 = bodyVertices[index++];
		vec3 v2 = bodyVertices[index++];

		// Get cross product of edges and normal vector.
		vec3 V1mV0 = v1 - v0;
		vec3 V2mV0 = v2 - v0;
		vec3 N = cross(V1mV0, V2mV0);

		// Compute integral terms.
		float tmp0, tmp1, tmp2;
		float f1x, f2x, f3x, g0x, g1x, g2x;
		tmp0 = v0.v[0] + v1.v[0];
		f1x = tmp0 + v2.v[0];
		tmp1 = v0.v[0] * v0.v[0];
		tmp2 = tmp1 + v1.v[0] * tmp0;
		f2x = tmp2 + v2.v[0] * f1x;
		f3x = v0.v[0] * tmp1 + v1.v[0] * tmp2 + v2.v[0] * f2x;
		g0x = f2x + v0.v[0] * (f1x + v0.v[0]);
		g1x = f2x + v1.v[0] * (f1x + v1.v[0]);
		g2x = f2x + v2.v[0] * (f1x + v2.v[0]);

		float f1y, f2y, f3y, g0y, g1y, g2y;
		tmp0 = v0.v[1] + v1.v[1];
		f1y = tmp0 + v2.v[1];
		tmp1 = v0.v[1] * v0.v[1];
		tmp2 = tmp1 + v1.v[1] * tmp0;
		f2y = tmp2 + v2.v[1] * f1y;
		f3y = v0.v[1] * tmp1 + v1.v[1] * tmp2 + v2.v[1] * f2y;
		g0y = f2y + v0.v[1] * (f1y + v0.v[1]);
		g1y = f2y + v1.v[1] * (f1y + v1.v[1]);
		g2y = f2y + v2.v[1] * (f1y + v2.v[1]);

		float f1z, f2z, f3z, g0z, g1z, g2z;
		tmp0 = v0.v[2] + v1.v[2];
		f1z = tmp0 + v2.v[2];
		tmp1 = v0.v[2] * v0.v[2];
		tmp2 = tmp1 + v1.v[2] * tmp0;
		f2z = tmp2 + v2.v[2] * f1z;
		f3z = v0.v[2] * tmp1 + v1.v[2] * tmp2 + v2.v[2] * f2z;
		g0z = f2z + v0.v[2] * (f1z + v0.v[2]);
		g1z = f2z + v1.v[2] * (f1z + v1.v[2]);
		g2z = f2z + v2.v[2] * (f1z + v2.v[2]);

		// Update integrals.
		integral[0] += N.v[0] * f1x;
		integral[1] += N.v[0] * f2x;
		integral[2] += N.v[1] * f2y;
		integral[3] += N.v[2] * f2z;
		integral[4] += N.v[0] * f3x;
		integral[5] += N.v[1] * f3y;
		integral[6] += N.v[2] * f3z;
		integral[7] += N.v[0] * (v0.v[1] * g0x + v1.v[1] * g1x + v2.v[1] * g2x);
		integral[8] += N.v[1] * (v0.v[2] * g0y + v1.v[2] * g1y + v2.v[2] * g2y);
		integral[9] += N.v[2] * (v0.v[0] * g0z + v1.v[0] * g1z + v2.v[0] * g2z);
	}

	integral[0] *= oneDiv6;
	integral[1] *= oneDiv24;
	integral[2] *= oneDiv24;
	integral[3] *= oneDiv24;
	integral[4] *= oneDiv60;
	integral[5] *= oneDiv60;
	integral[6] *= oneDiv60;
	integral[7] *= oneDiv120;
	integral[8] *= oneDiv120;
	integral[9] *= oneDiv120;

	// mass
	this->mass = integral[0];

	// center of mass
	this->bodyCOM = vec4(integral[1], integral[2], integral[3], 0.0f) / mass;

	// inertia relative to world origin
	this->Ibody.m[0] = integral[5] + integral[6];
	this->Ibody.m[1] = -integral[7];
	this->Ibody.m[2] = -integral[9];
	this->Ibody.m[3] = 0.0f;

	this->Ibody.m[4] = -integral[7];
	this->Ibody.m[5] = integral[4] + integral[6];
	this->Ibody.m[6] = -integral[8];
	this->Ibody.m[7] = 0.0f;

	this->Ibody.m[8] = -integral[9];
	this->Ibody.m[9] = -integral[8];
	this->Ibody.m[10] = integral[4] + integral[5];
	this->Ibody.m[11] = 0.0f;

	this->Ibody.m[12] = 0.0f;
	this->Ibody.m[13] = 0.0f;
	this->Ibody.m[14] = 0.0f;
	this->Ibody.m[15] = 1.0f;

	// inertia relative to center of mass
	if (bodyCoords)
	{
		this->Ibody.m[0] -= mass*(bodyCOM.v[1] * bodyCOM.v[1] + bodyCOM.v[2] * bodyCOM.v[2]);
		this->Ibody.m[1] += mass*bodyCOM.v[0] * bodyCOM.v[1];
		this->Ibody.m[2] += mass*bodyCOM.v[2] * bodyCOM.v[0];
		this->Ibody.m[3] = 0.0f;

		this->Ibody.m[4] += mass*bodyCOM.v[0] * bodyCOM.v[1];
		this->Ibody.m[5] -= mass*(bodyCOM.v[2] * bodyCOM.v[2] + bodyCOM.v[0] * bodyCOM.v[0]);
		this->Ibody.m[6] += mass*bodyCOM.v[1] * bodyCOM.v[2];
		this->Ibody.m[7] = 0.0f;

		this->Ibody.m[8] += mass*bodyCOM.v[2] * bodyCOM.v[0];
		this->Ibody.m[9] += mass*bodyCOM.v[1] * bodyCOM.v[2];
		this->Ibody.m[10] -= mass*(bodyCOM.v[0] * bodyCOM.v[0] + bodyCOM.v[1] * bodyCOM.v[1]);
		this->Ibody.m[11] = 0.0f;

		this->Ibody.m[12] = 0.0f;
		this->Ibody.m[13] = 0.0f;
		this->Ibody.m[14] = 0.0f;
		this->Ibody.m[15] = 1.0f;
	}
}

void RigidBody::drawMesh(mat4 view, mat4 projection, vec4 viewPosition)
{
	mat4 objectModel = scale(identity_mat4(), vec3(this->scaleFactor, this->scaleFactor, this->scaleFactor));
	objectModel = this->rotation * objectModel;
	objectModel = translate(objectModel, this->position);
	rigidBodyMesh.drawMesh(view, projection, objectModel, this->meshColour, viewPosition);
}

void RigidBody::drawBoundingSphere(mat4 view, mat4 projection)
{
	mat4 objectModel = scale(identity_mat4(), vec3(this->boundingSphereRadius, this->boundingSphereRadius, this->boundingSphereRadius));
	objectModel = this->rotation * objectModel;
	objectModel = translate(objectModel, this->position);
	boundingSphere.drawLine(view, projection, objectModel, this->boundingSphereColour);
}

void RigidBody::drawAABB(mat4 view, mat4 projection, GLuint* shaderID)
{
	GLfloat bounding_box_vertices[] = {
		xMax, yMax, zMax,
		xMin, yMax, zMax,
		xMin, yMin, zMax,
		xMax, yMin, zMax,

		xMax, yMax, zMax,
		xMax, yMax, zMin,
		xMax, yMin, zMin,
		xMax, yMin, zMax,

		xMin, yMin, zMax,
		xMin, yMin, zMin,
		xMax, yMin, zMin,
		xMax, yMax, zMin,

		xMin, yMax, zMin,
		xMin, yMin, zMin,
		xMin, yMax, zMin,
		xMin, yMax, zMax
	};

	Mesh bounding_box = Mesh(shaderID);
	bounding_box.generateObjectBufferMesh(bounding_box_vertices, 16);
	bounding_box.drawLine(view, projection, identity_mat4(), boundingBoxColour);

	bounding_box.dispose();
}
#pragma endregion

// Functions
float calculateImpulse(RigidBody &rbA, RigidBody &rbB, vec4 contactNormal, vec4 contactPoint, bool &collidingContact);
float calculatePlaneImpulse(RigidBody &rigidBody, vec4 contactPlaneNormal, vec4 contactPoint, bool &collidingContact);
void checkPlaneCollisions(RigidBody &rigidBody);

vec4 getTorque(vec4 force, vec4 position, vec4 point);
void computeForcesAndTorque(RigidBody &rigidBody);
void updateRigidBodies(GLuint mode, GLuint numRigidBodies, vector<RigidBody> &rigidbodies);

void checkBoundingSphereCollisions(GLuint numRigidBodies, vector<RigidBody> &rigidbodies);
bool isColliding(const RigidBody& bdi, const RigidBody& cdi);
void checkAABBCollisions(GLuint numRigidBodies, vector<RigidBody> &rigidbodies);

#pragma region COLLISION_RESPONSE
float calculateImpulse(RigidBody &rbA, RigidBody &rbB, vec4 contactNormal, vec4 contactPoint, bool &collidingContact)
{
	vec4 pA_dot = rbA.velocity + cross(rbA.angularVelocity, (contactPoint - rbA.worldCOM));
	vec4 pB_dot = rbB.velocity + cross(rbB.angularVelocity, (contactPoint - rbB.worldCOM));
	vec4 pA_pB = pA_dot - pB_dot;
	float v_relative = dot(vec3(contactNormal.v[0], contactNormal.v[1], contactNormal.v[2]), vec3(pA_pB.v[0], pA_pB.v[1], pA_pB.v[2]));

	float numerator = -(1 + restitution) * v_relative;

	float massA_inverse = 1 / rbA.mass;
	float massB_inverse = 1 / rbB.mass;

	vec4 rA = contactPoint - rbA.worldCOM;
	vec4 term3_a = cross(rA, contactNormal);
	vec4 term3_b = rbA.Iinv * term3_a;
	vec3 term3_c = cross(vec3(term3_b.v[0], term3_b.v[1], term3_b.v[2]), vec3(rA.v[0], rA.v[1], rA.v[2]));
	float term3 = dot(vec3(contactNormal.v[0], contactNormal.v[1], contactNormal.v[2]), term3_c);

	vec4 rB = contactPoint - rbB.worldCOM;
	vec4 term4_a = cross(rB, contactNormal);
	vec4 term4_b = rbB.Iinv * term4_a;
	vec3 term4_c = cross(vec3(term4_b.v[0], term4_b.v[1], term4_b.v[2]), vec3(rB.v[0], rB.v[1], rB.v[2]));
	float term4 = dot(vec3(contactNormal.v[0], contactNormal.v[1], contactNormal.v[2]), term4_c);

	float j_colliding = numerator / (massA_inverse + massB_inverse + term3 + term4);
	float j_resting = (-1 * v_relative) / (massA_inverse + massB_inverse + term3 + term4);

	if (v_relative >= 0.001f) //Moving away
	{
		collidingContact = false;
		return 0.0f;
	}
	else if (v_relative >= -0.001f) // Resting Contact
	{
		collidingContact = false;
		return j_resting;
	}
	else
	{
		// Colliding Contact
		collidingContact = true;
		return j_colliding;
	}
}

float calculatePlaneImpulse(RigidBody &rigidBody, vec4 contactPlaneNormal, vec4 contactPoint, bool &collidingContact)
{
	vec4 pA_dot = rigidBody.velocity + cross(rigidBody.angularVelocity, (contactPoint - rigidBody.worldCOM));
	float v_relative = dot(vec3(contactPlaneNormal.v[0], contactPlaneNormal.v[1], contactPlaneNormal.v[2]), vec3(pA_dot.v[0], pA_dot.v[1], pA_dot.v[2]));

	float modifier = -(1.0f + restitution);

	float numerator = modifier * v_relative;

	float mass_inverse = 1 / rigidBody.mass;

	vec4 rA = contactPoint - rigidBody.worldCOM;
	vec4 term3_a = cross(rA, contactPlaneNormal);
	vec4 term3_b = rigidBody.Iinv * term3_a;
	vec3 term3_c = cross(vec3(term3_b.v[0], term3_b.v[1], term3_b.v[2]), vec3(rA.v[0], rA.v[1], rA.v[2]));
	float term3 = dot(vec3(contactPlaneNormal.v[0], contactPlaneNormal.v[1], contactPlaneNormal.v[2]), term3_c);

	float j_colliding = numerator / (mass_inverse + term3);
	float j_resting = (-1 * v_relative) / (mass_inverse + term3);

	if (v_relative >= 0.001f) //Moving away
	{
		collidingContact = false;
		return 0.0f;
	}
	else if (v_relative >= -0.001f) // Resting Contact
	{
		collidingContact = false;
		return j_colliding;
	}
	else
	{
		// Colliding Contact
		collidingContact = true;
		return j_colliding;
	}

}

void checkPlaneCollisions(RigidBody &rigidBody)
{
	// Add all of the plane normals
	vector<vec4> plane_normals;
	plane_normals.push_back(yAxis);
	plane_normals.push_back(yAxis * -1);
	plane_normals.push_back(xAxis);
	plane_normals.push_back(xAxis * -1);
	plane_normals.push_back(zAxis);
	plane_normals.push_back(zAxis * -1);

	for (int i = 0; i < (int)plane_normals.size(); i++)
	{
		if (pointToPlane(rigidBody.worldCentroid, plane_normals[i], (plane_normals[i] * -1)) <= rigidBody.boundingSphereRadius)
		{
			for (GLuint j = 0; j < rigidBody.numPoints; j++)
			{
				float distanceToPlane = pointToPlane(rigidBody.worldVertices[j], plane_normals[i], (plane_normals[i] * -1));
				if (distanceToPlane < 0.001f)
				{
					bool collidingContact = false;
					float impulseMagnitude = calculatePlaneImpulse(rigidBody, plane_normals[i], rigidBody.worldVertices[j], collidingContact);
					//printf("Impulse magnitude: %f \n", impulseMagnitude);
					vec4 impulseForce = plane_normals[i] * impulseMagnitude;
					
					rigidBody.force = impulseForce;
					rigidBody.linearMomentum += rigidBody.force;
					rigidBody.velocity = rigidBody.linearMomentum / rigidBody.mass;
					vec4 normalVelocity = plane_normals[i] * dot(vec3(plane_normals[i].v[0], plane_normals[i].v[1], plane_normals[i].v[2]), vec3(rigidBody.velocity.v[0], rigidBody.velocity.v[1], rigidBody.velocity.v[2]));
					if (vec4Magnitude(normalVelocity) < 0.1f)
						rigidBody.velocity = rigidBody.velocity - normalVelocity;
					
					//if (vec4Magnitude(rigidBody.velocity) <= 0.1f)
					//	rigidBody.velocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);
					
					//if (collidingContact)
					//{
						rigidBody.torque = getTorque(impulseForce, rigidBody.worldCOM, rigidBody.worldVertices[j]);
						rigidBody.angularMomentum += rigidBody.torque;	
						rigidBody.angularVelocity = rigidBody.Iinv * rigidBody.angularMomentum;
						//if (vec4Magnitude(rigidBody.angularVelocity) <= 0.1f)
						//	rigidBody.angularVelocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);
					//}		
				}
			}

			// Previous Implementation
			/*rigidBody.velocity = rigidBody.linearMomentum / rigidBody.mass;
			if (vec4Magnitude(rigidBody.velocity) <= 0.1f)
				rigidBody.velocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);

			rigidBody.angularVelocity = rigidBody.Iinv * rigidBody.angularMomentum;
			if (vec4Magnitude(rigidBody.angularVelocity) <= 0.1f)
				rigidBody.angularVelocity = vec4(0.0f, 0.0f, 0.0f, 0.0f);*/
		}
	}
}
#pragma endregion

#pragma region RIGIDBODY_UPDATE
vec4 getTorque(vec4 force, vec4 position, vec4 point)
{
	vec4 pointToCOM = point - position;
	return cross(pointToCOM, force); 
}

void computeForcesAndTorque(RigidBody &rigidBody)
{
	// Clear Forces
	rigidBody.force = vec4(0.0f, 0.0f, 0.0f, 0.0f);
	rigidBody.torque = vec4(0.0f, 0.0f, 0.0f, 0.0f);

	// Apply Gravity
	if(useGravity)
		rigidBody.force += gravity * 1.0f * rigidBody.mass;

	if (useForce == FORCE1)
	{
		rigidBody.force -= gravity * 0.01f * rigidBody.mass;
	}
	else if (useForce == FORCE2)
	{
		rigidBody.force += vec4(0.001f, 0.001f, 0.001f, 0.0f) * rigidBody.mass;
	}
	else if (useForce == FORCE3)
	{
		rigidBody.force += vec4(-0.001f, 0.001f, 0.001f, 0.0f) * rigidBody.mass;
	}
	else if (useForce == FORCE4)
	{
		rigidBody.force += vec4(0.001f, 0.0f, 0.001f, 0.0f) * rigidBody.mass;
	}
}


void updateRigidBodies(GLuint mode, GLuint numRigidBodies, vector<RigidBody> &rigidbodies)
{
	for (GLuint i = 0; i < numRigidBodies; i++)
	{
		RigidBody &rigidBody = rigidbodies[i];

		computeForcesAndTorque(rigidBody);

		rigidBody.position += rigidBody.velocity * deltaTime;

		versor omega;
		omega.q[0] = 0.0f;
		omega.q[1] = rigidBody.angularVelocity.v[0];
		omega.q[2] = rigidBody.angularVelocity.v[1];
		omega.q[3] = rigidBody.angularVelocity.v[2];

		versor angularVelocityQuat;
		float avMag = quatMagnitude(omega);
		angularVelocityQuat.q[0] = cos((avMag * deltaTime) / 2);
		if (avMag > 0)
		{
			angularVelocityQuat.q[1] = (rigidBody.angularVelocity.v[0] / avMag) * sin((avMag * deltaTime) / 2);
			angularVelocityQuat.q[2] = (rigidBody.angularVelocity.v[1] / avMag) * sin((avMag * deltaTime) / 2);
			angularVelocityQuat.q[3] = (rigidBody.angularVelocity.v[2] / avMag) * sin((avMag * deltaTime) / 2);
		}
		else
		{
			angularVelocityQuat.q[1] = 0.0f;
			angularVelocityQuat.q[2] = 0.0f;
			angularVelocityQuat.q[3] = 0.0f;
		}

		multiplyQuat(rigidBody.orientation, angularVelocityQuat, rigidBody.orientation);

		rigidBody.linearMomentum += rigidBody.force * deltaTime;
		rigidBody.angularMomentum += rigidBody.torque * deltaTime;

		rigidBody.velocity = rigidBody.linearMomentum / rigidBody.mass;
		rigidBody.rotation = quat_to_mat4(normalise(rigidBody.orientation));
		rigidBody.Iinv = rigidBody.rotation * rigidBody.IbodyInv * transpose(rigidBody.rotation);
		rigidBody.angularVelocity = rigidBody.Iinv * rigidBody.angularMomentum;

		// Update all world points
		for (GLuint i = 0; i < rigidBody.numPoints; i++)
		{
			rigidBody.worldVertices[i] = (rigidBody.rotation * rigidBody.initialWorldVertices[i]) + rigidBody.position;
		}

		rigidBody.worldCOM = (rigidBody.rotation * rigidBody.bodyCOM) + rigidBody.position;
		rigidBody.worldCentroid = (rigidBody.rotation * rigidBody.bodyCentroid) + rigidBody.position;

		checkPlaneCollisions(rigidBody);

		if (mode == AABB)
		{
			rigidBody.xMin = rigidBody.worldVertices[0].v[0];
			rigidBody.xMax = rigidBody.worldVertices[0].v[0];
			rigidBody.yMin = rigidBody.worldVertices[0].v[1];
			rigidBody.yMax = rigidBody.worldVertices[0].v[1];
			rigidBody.zMin = rigidBody.worldVertices[0].v[2];
			rigidBody.zMax = rigidBody.worldVertices[0].v[2];

			for (GLuint i = 1; i < rigidBody.numPoints; i++)
			{
				vec4 vertex = rigidBody.worldVertices[i];

				if (vertex.v[0] < rigidBody.xMin)
					rigidBody.xMin = vertex.v[0];
				else if (vertex.v[0] > rigidBody.xMax)
					rigidBody.xMax = vertex.v[0];

				if (vertex.v[1] < rigidBody.yMin)
					rigidBody.yMin = vertex.v[1];
				else if (vertex.v[1] > rigidBody.yMax)
					rigidBody.yMax = vertex.v[1];

				if (vertex.v[2] < rigidBody.zMin)
					rigidBody.zMin = vertex.v[2];
				else if (vertex.v[2] > rigidBody.zMax)
					rigidBody.zMax = vertex.v[2];
			}
		}

		// Reset the colliding with counter
		rigidBody.collidingWith = 0;
	}
}
#pragma endregion

#pragma region BROAD_PHASE_COLLISION_DETECTION
void checkBoundingSphereCollisions(GLuint numRigidBodies, vector<RigidBody> &rigidbodies)
{
	for (GLuint i = 0; i < numRigidBodies; i++)
	{
		RigidBody &rb1 = rigidbodies[i];

		for (GLuint j = i + 1; j < numRigidBodies; j++)
		{
			RigidBody &rb2 = rigidbodies[j];

			if (getDistance(rb1.position, rb2.position) <= (rb1.boundingSphereRadius + rb2.boundingSphereRadius))
			{
				rb1.collidingWith++;
				rb2.collidingWith++;
				rb1.boundingSphereColour = red;
				rb2.boundingSphereColour = red;
			}
		}

		if (rb1.collidingWith == 0)
		{
			rb1.boundingSphereColour = green;
		}
	}
}

bool isColliding(const RigidBody& bdi, const RigidBody& cdi)
{
	if ((bdi.xMax < cdi.xMin || bdi.xMin > cdi.xMax)) return false;
	if ((bdi.yMax < cdi.yMin || bdi.yMin > cdi.yMax)) return false;
	if ((bdi.zMax < cdi.zMin || bdi.zMin > cdi.zMax)) return false;

	return true;
}

void checkAABBCollisions(GLuint numRigidBodies, vector<RigidBody> &rigidbodies)
{
	vector<bool> m_collision(numRigidBodies, false);
	for (GLuint i = 0; i < numRigidBodies; i++)
	{
		for (GLuint j = i + 1; j < numRigidBodies; j++)
		{
			if (isColliding(rigidbodies[i], rigidbodies[j]))
			{
				m_collision[i] = true;
				m_collision[j] = true;
			}
		}
	}

	for (GLuint k = 0; k < numRigidBodies; k++)
	{
		rigidbodies[k].collisionAABB = m_collision[k];
		rigidbodies[k].boundingBoxColour = m_collision[k] ? red : green;
	}
}
#pragma endregion