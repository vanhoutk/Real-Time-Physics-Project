#pragma once

#include "Antons_maths_funcs.h"

#define INFINITY	1.0e+20
#define EPSILON		0.0001
#define CONE_EPS	1.0e-7

/*
 * Flaring Angles
 * 
 * These are the angles by which to flare the various planes defining
 * Voronoi cones.  Planes of vertex Voronoi regions are flared in,
 * planes of face Voronoi regions are flared out, and planes of edge
 * Voronoi regions are flared out if they correspond to a vertex on the
 * edge's coboundary, or in if they correspond to a face on the coboundary.
 * Stated simply, the highest dimensional feature is given priority.

 * All these constant are in radians and should be nonnegative.  Also, the
 * following inequalities must hold:

 * FACE_FLARE > EDGE_FACE_FLARE
 * FACE_FLARE > VERTEX_FLARE
 * EDGE_VERTEX_FLARE > VERTEX_FLARE
*/

#define VERTEX_FLARE		0.5e-5
#define EDGE_VERTEX_FLARE	1.0e-5
#define EDGE_FACE_FLARE		1.0e-5
#define FACE_FLARE			2.0e-5

/*
=============================================================================

data structures

1. VERTEX, EDGE, and FACE

These are the structures for vertices, edges, and faces.  The fields which
appear in these structures are listed below.

tag:  This is an integer field which appears first in all three structures.
It takes a value V, E, or F, and is provided so we can identify the
type of an arbitrary feature (e.g. one pointed to by a pointer to void)

name:  This string field is present in the VERTEX and FACE structures.
Edge names are derived from the names of the endpoint vertices.

coords & xcoords:  coords are the coordinates of a VERTEX.  When we
transform the VERTEX to a different frame, xcoords holds the transformed
coordinates.

u & xu:  u is a unit vector pointing from endpoint v1 to endpoint v2 of
an EDGE.  xu is the transformed version when we transform the edge to
a different frame.

len:  The length of an EDGE.

plane:  This field is present in the FACE structure.  It defines the
four coefficients of the plane of the face, i.e. the plane of the
face is given by plane[0]*x + plane[1]*y + plane[2]*z + plane[3] = 0.
Furthermore, the first 3 components are a unit vector, which form
the outward normal of the face.

verts:  This is a list of vertices (in CCW order) around a FACE.

edges:  This is a list of edges incident to a VERTEX, or the boundary
edges (in CCW order) of a FACE.

v1 & v2:  The startpoint end endpoint vertices of an EDGE.

fl & fr:  The left and right faces of an EDGE, respectively.  The
"left" face is the one on the left hand side of an observer standing
on v1, looking toward v2, and standing on the outside of the polyhedron.

cone:  This field is present in all three structures, and points to
the list of planes defining the features voronoi region.  For a
FACE, only the side planes are present in this list.


2.  POLYHEDRON

This is the basic data structure for a polyhedron.  It contains the
following fields.

name:  A string; the name of the polyhedron, if any.

verts:  The list of vertices of the polyhedron.

edges:  The list of edges of the polyhedron.

faces:  What a surprise -- the list of faces of the polyhedron.

pose:  The pose field is a transformation matrix.  The coordinates of
a polyhedron's vertices, and hence the location of its edges and faces
are all defined relative to a body frame.  All of the planes defining
the voronoi region are also relative to this frame.  However, a
polyhedron is allowed to move through space, and the pose matrix
expresses the pose of the body frame relative to some fixed global
frame.  Thus, the pose matrix Px of polyhedron x transforms body
coordinates to "world" coordinates.  Now let x and y be two polyhedra,
with pose matrices Px and Py.  Suppose we want to compute Txy, a
transformation matrix which transforms coordinates from x's body frame
to y's body frame, and Tyx, it's inverse.  It's easy to show that
Txy = Inv(Py) * Px, and Tyx = Inv(Px) * Py.


3. FNODE

FNODE stands for feature node.  FNODEs are used for building lists of
features, such as the lists of features in the POLYHEDRON structure, and
also the list of edges in the VERTEX structure, and the lists of vertices
and edges in the FACE structure.  An FNODE has two fields.

f:  f points to an actual feature.  f is a union, so it may be treated
as a pointer to a VERTEX, EDGE, or FACE, depending on what type of
feature list we're dealing with.  There is also an "any" field in the
union, for when we don't care what type of feature f points to.

next:  This points to the next FNODE in the list, or is NULL if this is
the last one.


4. PNODE

PNODE stands for plane node.  PNODEs are used for building lists of
planes which define the voronoi regions for the features of a polyhedron.
These lists are pointed to by the cone fields of the VERTEX, EDGE, and
FACE structure.  PNODEs contain the following three fields.

plane:  These are the four coefficients of the actual plane, i.e.
the plane is given by plane[0]*x + plane[1]*y + plane[2]*z + plane[3] = 0.
Furthermore, the first three components of plane form a unit vector
pointing toward the "good" side of the plane.  If a point lies within
a voronoi region, it lies on the good side of all planes defining that
region.  If it lies on the bad sign of any of the planes, it is not in
the voronoi region.

nbr:  This is points to the neighboring feature which gives rise to
the particular plane specified by the plane field.  When a point
lies on the "bad" side of a particular plane, the nbr field points
to the feature we should walk to when searching for the closest
feature pair.

next:  This points to the next PNODE in the list, or is NULL if this is
the last one.


=============================================================================
*/

struct vertex 
{
	int tag;			// Used for identifying feature type (V, E, F)
	char name[20];
	vec3 coords;		// Coordinates
	vec3 xcoords;		// Tranformed Coordinates in word space
	featureNode *edges;	// List of edges incident to a vertex
	planeNode *cone;	// Points to the list of planes defining the feature's voronoi region
};
typedef struct vertex VERTEX;

struct edge 
{
	int tag;		// Used for identifying feature type (V, E, F)
	vertex *v1, *v2;// Startpoint and endpoint vertices of an edge
	face *fl, *fr;  // Face to the left and face to the right of the edge
	vec3 u;			// unit vector from v1 -> v2 
	vec3 xu;		// transformed unit vector from v1 -> v2 
	float len;		// length of edge 
	planeNode *cone;// Points to the list of planes defining the feature's voronoi region
};
typedef struct edge EDGE;

struct face 
{
	int tag;				// Used for identifying feature type (V, E, F)
	char name[20];
	featureNode *vertices;	// List of vertices in CCW order around a face
	featureNode *edges;		// List of boundary edges in CCW order of a face
	vec4 plane;				// It defines the four coefficients of the plane of the face, i.e.the plane of the
							// face is given by plane[0] * x + plane[1] * y + plane[2] * z + plane[3] = 0.
							// Furthermore, the first 3 components are a unit vector, which form
							// the outward normal of the face.
	planeNode *cone;		// Points to the list of planes defining the feature's voronoi region
							// For a face, only the side planes are present in this list.
};
typedef struct face FACE;

struct polyhedron 
{
	char name[20];		// Name of the polyhedron, if any
	mat4 pose;			// Transformation of the polyhedron
	featureNode *verts; // List of vertices of the polyhedron
	featureNode *edges; // List of edges of the polyhedron
	featureNode *faces; // List of faces of the polyhedron
};
typedef struct polyhedron POLYHEDRON;

struct featureNode 
{
	union 
	{
		vertex *v;
		edge *e;
		face *f;
		void *any;
	} f;					  // Points to an actual feature of any of the above types
	struct featureNode *next; // Next feature node in the list, or null if the last one 
};
typedef struct featureNode FNODE;

struct planeNode 
{
	float plane[4];
	// plane[0] * x + plane[1] * y + plane[2] * z + plane[3] >= 0 
	void *nbr; // point to feature to move to if this test fails 
	struct planeNode *next;  // if there are more planes in this cone
};
typedef struct planeNode PNODE;

/*
 * Feature Types
 *
 * These are the values stored in the tag fields of the vertex, edge, and
 * face structures.  The featTag macro returns the type of an arbitrary feature.
 */

#define V 1
#define E 2
#define F 3
#define featureTag(featurePtr) (*((int *) (featurePtr)))

/*
 * Memory Allocation Macros
 */

#define allocVertex (vertex *) malloc(sizeof(vertex))
#define allocEdge (edge *) malloc(sizeof(edge))
#define allocFace (face *) malloc(sizeof(face))
#define allocPolyhedron (polyhedron *) malloc(sizeof(polyhedron))
#define allocFnode (featureNode *) malloc(sizeof(featureNode))
#define allocPnode (planeNode *) malloc(sizeof(planeNode))

#pragma region NOT_USED
/*
 * Transformation Macros
 *
 * These transform a vertex or an edge from one frame to another using the
 * transformation matrix T.  Transforming an vertex involves transforming
 * a single set of coordinates.  Transforming an edge involves transforming
 * the coordinates of both its endpoints, as well as its direction vector
 */

// #define xformVert(T, v) xformPoint(T, ((VERTEX *)(v))->coords, ((VERTEX *)(v))->xcoords)
// #define xformEdge(T, e) {xformPoint(T, ((EDGE *)(e))->v1->coords, ((EDGE *)(e))->v1->xcoords); xformPoint(T, ((EDGE *)(e))->v2->coords, ((EDGE *)(e))->v2->xcoords); xformVect(T, ((EDGE *)(e))->u, ((EDGE *)(e))->xu);}

/*
 * Set the pose matrix of polyhedron p
 */

// #define setPose(p, T) mat4copy(T, (p)->pose)
#pragma endregion

/*
 * Miscellaneous macros
 */

#define SWAP(x,y,type) {type tmp; tmp = x; x = y; y = tmp;}
#define SQR(x) ((x) * (x))
#define flush fflush(stdout)

/*
 * Function declarations
 */

int getWord();
int random2();
void addFeature();
void reverseFlist();
void featureName();
VERTEX *findVertex();
EDGE *findEdge();
FACE *findFace();
int numFeats();
void *nthFeat();
void *randFeat();
float polyRad();
EDGE *newEdge();
int loadPolyhedronLibrary();
POLYHEDRON *createPolyhedron();
void dumpPolyhedron();
void addPlane();
void flipPlane();
void tweakPlaneNorm();
void computeVertCone();
void computeFaceCone();
void buildCones();
void dumpCones();
int vertConeChk();
int edgeConeChk();
int faceConeChk();
float Dvv();
float Dev();
float Dfv();
float Dve();
float Dee();
float Dfe();
float Dff();
void *closestToVert();
void *closestToEdge();
void *closestToFace();
void edgeCPs();
void *closestEdgeOrVertOnFace();
void closestEdges();
void *closestToFacePlane();
int polygonCut();
int faceOverlap();
float vertex_vertex();
float vertex_edge();
float edge_edge();
float vertex_face();
float edge_face();
float face_face();
float closestFeatures();
float closestFeaturesInit();