#pragma once


#include "Antons_maths_funcs.h"
//#define INFINITY	(float)1.0e+20
#define EPSILON		(float)0.0001
#define CONE_EPS	(float)1.0e-7


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

#define VERTEX_FLARE		(float)0.5e-5
#define EDGE_VERTEX_FLARE	(float)1.0e-5
#define EDGE_FACE_FLARE		(float)1.0e-5
#define FACE_FLARE			(float)2.0e-5

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
	struct featureNode *edges;	// List of edges incident to a vertex
	struct planeNode *cone;	// Points to the list of planes defining the feature's voronoi region
};
typedef struct vertex VERTEX;

struct edge 
{
	int tag;		// Used for identifying feature type (V, E, F)
	struct vertex *v1, *v2;// Startpoint and endpoint vertices of an edge
	struct face *fl, *fr;  // Face to the left and face to the right of the edge
	vec3 u;			// unit vector from v1 -> v2 
	vec3 xu;		// transformed unit vector from v1 -> v2 
	float len;		// length of edge 
	struct planeNode *cone;// Points to the list of planes defining the feature's voronoi region
};
typedef struct edge EDGE;

struct face 
{
	int tag;			// Used for identifying feature type (V, E, F)
	char name[20];
	struct featureNode *verts;	// List of vertices in CCW order around a face
	struct featureNode *edges;	// List of boundary edges in CCW order of a face
	struct vec4 plane;			// It defines the four coefficients of the plane of the face, i.e.the plane of the
						// face is given by plane[0] * x + plane[1] * y + plane[2] * z + plane[3] = 0.
						// Furthermore, the first 3 components are a unit vector, which form
						// the outward normal of the face.
	struct planeNode *cone;	// Points to the list of planes defining the feature's voronoi region
						// For a face, only the side planes are present in this list.
};
typedef struct face FACE;

struct polyhedron 
{
	char name[20];		// Name of the polyhedron, if any
	struct mat4 pose;			// Transformation of the polyhedron
	struct featureNode *verts; // List of vertices of the polyhedron
	struct featureNode *edges; // List of edges of the polyhedron
	struct featureNode *faces; // List of faces of the polyhedron
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
	struct vec4 plane;
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

#define setPose(p, T) (p)->pose = T


#pragma region XFORM_FUNCTIONS
void xformPoint(mat4 M, vec3 p, vec3 &result)
{
	for (int i = 0; i < 3; i++)
		result.v[i] = M.m[i] * p.v[0] + M.m[i + 4] * p.v[1] + M.m[i + 8] * p.v[2] + M.m[i + 12];

	mat3 rotation;
	rotation.m[0] = M.m[0];
	rotation.m[1] = M.m[1];
	rotation.m[2] = M.m[2];
	rotation.m[3] = M.m[4];
	rotation.m[4] = M.m[5];
	rotation.m[5] = M.m[6];
	rotation.m[6] = M.m[8];
	rotation.m[7] = M.m[9];
	rotation.m[8] = M.m[10];

	vec3 position = vec3(M.m[12], M.m[13], M.m[14]);

	vec3 result2 = (rotation * p) + position;
}

void xformVector(mat4 M, vec3 p, vec3 &result)
{
	for (int i = 0; i < 3; i++)
		result.v[i] = M.m[i] * p.v[0] + M.m[i + 4] * p.v[1] + M.m[i + 8] * p.v[2];

	mat3 rotation;
	rotation.m[0] = M.m[0];
	rotation.m[1] = M.m[1];
	rotation.m[2] = M.m[2];
	rotation.m[3] = M.m[4];
	rotation.m[4] = M.m[5];
	rotation.m[5] = M.m[6];
	rotation.m[6] = M.m[8];
	rotation.m[7] = M.m[9];
	rotation.m[8] = M.m[10];

	vec3 result2 = (rotation * p);
}

void xformVertex(mat4 M, vertex *v)
{
	xformPoint(M, v->coords, v->xcoords);
}

void xformEdge(mat4 M, edge *e)
{
	xformPoint(M, e->v1->coords, e->v1->xcoords);
	xformPoint(M, e->v2->coords, e->v2->xcoords);
	xformVector(M, e->u, e->xu);
}

/*
 * Transformation matrix inversion:  Inverse(M) => inv
 * N.B. M and inv should not point to the same matrix.  We assume M is a
 * transformation matrix; this will not properly invert arbitrary 4x4 matrices.
 */

void matInvXform(mat4 M, mat4 &inv)
{
	// We invert the rotation part by transposing it 
	inv.m[0] = M.m[0];
	inv.m[1] = M.m[4];
	inv.m[2] = M.m[8];
	inv.m[4] = M.m[1];
	inv.m[5] = M.m[5];
	inv.m[6] = M.m[9];
	inv.m[8] = M.m[2];
	inv.m[9] = M.m[6];
	inv.m[10] = M.m[10];

	// The new displacement vector is given by:  d' = -(R^-1) * d 
	inv.m[12] = -inv.m[0] * M.m[12] - inv.m[4] * M.m[13] - inv.m[8] * M.m[14];
	inv.m[13] = -inv.m[1] * M.m[12] - inv.m[5] * M.m[13] - inv.m[9] * M.m[14];
	inv.m[14] = -inv.m[2] * M.m[12] - inv.m[6] * M.m[13] - inv.m[10] * M.m[14];

	// The rest stays the same
	inv.m[3] = inv.m[7] = inv.m[11] = 0.0;
	inv.m[15] = 1.0;
}

/*
 * Transformation matrix multiply:  a * b => c
 * N.B. This routine is much faster than the general 4 x 4 matrix
 * multiply above, but only works properly if a and b are SE(3) transformation
 * matrices.  c should not point to the same matrix as a or b!
 */
void matMultXform(mat4 a, mat4 b, mat4 &c)
{
	// Rc = Ra Rb
	c.m[0] = a.m[0] * b.m[0] + a.m[4] * b.m[1] + a.m[8] * b.m[2];
	c.m[4] = a.m[0] * b.m[4] + a.m[4] * b.m[5] + a.m[8] * b.m[6];
	c.m[8] = a.m[0] * b.m[8] + a.m[4] * b.m[9] + a.m[8] * b.m[10];
	c.m[1] = a.m[1] * b.m[0] + a.m[5] * b.m[1] + a.m[9] * b.m[2];
	c.m[5] = a.m[1] * b.m[4] + a.m[5] * b.m[5] + a.m[9] * b.m[6];
	c.m[9] = a.m[1] * b.m[8] + a.m[5] * b.m[9] + a.m[9] * b.m[10];
	c.m[2] = a.m[2] * b.m[0] + a.m[6] * b.m[1] + a.m[10] * b.m[2];
	c.m[6] = a.m[2] * b.m[4] + a.m[6] * b.m[5] + a.m[10] * b.m[6];
	c.m[10] = a.m[2] * b.m[8] + a.m[6] * b.m[9] + a.m[10] * b.m[10];

	// Vc = Ra Vb + Va
	c.m[12] = a.m[0] * b.m[12] + a.m[4] * b.m[13] + a.m[8] * b.m[14] + a.m[12];
	c.m[13] = a.m[1] * b.m[12] + a.m[5] * b.m[13] + a.m[9] * b.m[14] + a.m[13];
	c.m[14] = a.m[2] * b.m[12] + a.m[6] * b.m[13] + a.m[10] * b.m[14] + a.m[14];

	// The rest
	c.m[3] = c.m[7] = c.m[11] = 0.0;
	c.m[15] = 1.0;
}
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

int getWord(FILE *fp, char *s);
int randomNumber(int n);

void addFeature(void *feature, featureNode **list);
void reverseFlist(featureNode **list);
void featureName(void *feature, char *name);
vertex *findVertex(polyhedron *p, char *name);
vertex *findVertex(polyhedron *p, vec3 coords);
edge *findEdge(polyhedron *p, char *name1, char *name2);
face *findFace(polyhedron *p, char *name);
int numFeatures(polyhedron *p);
void *nthFeature(polyhedron *p, int n);
void *randFeat(polyhedron *p);
float polyhedronRadius(char *name);
edge *newEdge(vertex *v1, vertex *v2);

//int loadPolyhedronLibrary();
//polyhedron *createPolyhedron();
polyhedron createPolyhedron(char *name, vector<vec4> bodyVertices, mat4 transformation);
void dumpPolyhedron(polyhedron *p);

void addPlane(planeNode *pn, planeNode **cone);
void flipPlane(vec4 src, vec4 &dest);
void tweakPlaneNormal(vec3 u, float epsilon, vec3 nOrig, vec4 &nTweak);
void computeVertexCone(vertex *v);
void computeFaceCone(face *f);
void buildCones(polyhedron *p);
void dumpCones(polyhedron *p);

int vertexConeCheck(vertex **v, vec3 point, int update);
int edgeConeCheck(edge **e, vec3 point, int update);
int faceConeCheck(face **f, vec3 point, int update);

float Dvv(vertex *v1, vertex *v2);
float Dev(edge *e, vertex *v);
float Dfv(face *f, vertex *v);
float Dve(vertex *v, edge *e);
float Dee(edge *e1, edge *e2);
float Dfe(face *f, edge *e);
float Dff(face *f1, face *f2, mat4 T12, mat4 T21);

void *closestToVertex(vertex *v, polyhedron *p, mat4 Tvp);
void *closestToEdge(edge *e, polyhedron *p, mat4 Tep);
void *closestToFace(face *f, polyhedron *p, mat4 Tfp, mat4 Tpf);
void edgeCPs(edge **e1, edge **e2, vec3 &cp1, vec3 &cp2);
void *closestEdgeOrVertexOnFace(face *f, edge *e);
void closestEdges(face *f1, face *f2, mat4 T12, edge **closest1, edge **closest2);
void *closestToFacePlane(face *f1, face *f2, mat4 T12, mat4 T21);
int polygonCut(face *f, edge *e, float *min, float *max);
int faceOverlap(face *f1, face *f2, mat4 T12, mat4 T21);

float vertex_vertex(vertex **v1, vertex **v2, mat4 T12, mat4 T21);
float vertex_edge(vertex **v, edge **e, mat4 T12, mat4 T21);
float edge_edge(edge **e1, edge **e2, mat4 T12, mat4 T21);
float vertex_face(vertex **v, face **f, mat4 T12, mat4 T21, polyhedron *facePoly);
float edge_face(edge **e, face **f, mat4 T12, mat4 T21, polyhedron *facePoly);
float face_face(face **f1, face **f2, mat4 T12, mat4 T21, polyhedron *face1poly, polyhedron *face2poly);
float closestFeatures(polyhedron *poly1, void **feat1, polyhedron *poly2, void **feat2);
float closestFeaturesInit(polyhedron *poly1, void **feat1, polyhedron *poly2, void **feat2);




//#include <string.h>
//#include <stdlib.h>
////#include "Antons_maths_funcs.h"
//#include "ClosestFeatures.h"



/* CLoset features */

/*
* Global Variables
*/

polyhedron polyhedronLibrary[20];
int polyhedronLibraryCount;

/*
=============================================================================

Miscellaneous routines

These are general functions which are not specifically related to the
distance algorithm, but are called by some of the routines.

=============================================================================
*/


/*
* Get next word on the line; skip over spaces and tabs.
* return true if we read a word, false if we hit \n first.
*/

int getWord(FILE *fp, char *s)
{
	char c;

	do fscanf_s(fp, "%c", &c, 1); while ((c == ' ') || c == '\t');
	if (c == '\n') return 0;
	ungetc(c, fp);
	fscanf_s(fp, "%s", s, sizeof(s));
	return 1;
}


/*
* return random # from {0, 1, ..., n-1}
*/
int randomNumber(int n)
{
	float x;
	x = ((float)rand()) / 0x7FFFFFFF;  /* x = random # in interval [0, 1) */
	return (int)floor(x * n);
}

/*
=============================================================================

Polyhedra and feature housekeeping

These routines are for accessing the polyhedron library, and performing
basic operations on polyhedra and their features.

=============================================================================
*/

/*
* Add a feature to a feature list.  This is a generic routine callable for
* any type of feature and feature list.
*/
void addFeature(void *feature, featureNode **list)
{
	featureNode *feature_node;

	feature_node = allocFnode;
	feature_node->f.any = feature;
	feature_node->next = *list;
	*list = feature_node;
}

/*
* Reverse feature list

* Because of the way feature lists are constructed (in a stack-like manner),
* it is sometimes necessary to reverse the order of the list.  This is
* necessary, for example, to make the vertices in a face's vertex list
* appear in the list in CCW order.  This routine reverse any feature list.
*/
void reverseFlist(featureNode **list)
{
	featureNode *last, *feature_node, *next;

	for (last = NULL, feature_node = *list; feature_node; feature_node = next)
	{
		next = feature_node->next;
		feature_node->next = last;
		last = feature_node;
	}

	if (last) *list = last;
}

/*
* Get name of an arbitrary feature f, and return in name.  Vertices and faces
* have explicit names; the name of an edge is derived from the names of
* its endpoint vertices (in alphabetical order), separated by a "." .
*/
void featureName(void *feature, char *name)
{
	int type;

	type = featureTag(feature);
	if (type == V) strcpy_s(name, sizeof(name), ((vertex *)feature)->name);
	else if (type == F) strcpy_s(name, sizeof(name), ((face *)feature)->name);
	else sprintf_s(name, 2 * sizeof(name), "%s.%s", ((edge *)feature)->v1->name, ((edge *)feature)->v2->name);
}

/*
* Return pointer to a vertex stucture, given the vertex's name and the
* polyhedron where it lives.
*/
vertex *findVertex(polyhedron *p, char *name)
{
	featureNode *feature_node;

	for (feature_node = p->verts; feature_node && strcmp(name, feature_node->f.v->name); feature_node = feature_node->next);
	return (feature_node) ? feature_node->f.v : NULL;
}

vertex *findVertex(polyhedron *p, vec3 coords)
{
	featureNode *feature_node;

	for (feature_node = p->verts; feature_node && feature_node->f.v->coords != coords; feature_node = feature_node->next);
	return (feature_node) ? feature_node->f.v : NULL;
}

/*
* Return pointer to an edge stucture, given the names of the edges endpoint
* vertices, and the polyhedron where it lives.  The names of the vertices
* may be passed in either order.
*/
edge *findEdge(polyhedron *p, char *name1, char *name2)
{
	featureNode *feature_node;

	if (strcmp(name1, name2) > 0) SWAP(name1, name2, char *);
	for (feature_node = p->edges; feature_node &&
		(strcmp(name1, feature_node->f.e->v1->name) || strcmp(name2, feature_node->f.e->v2->name));
		feature_node = feature_node->next);
	return (feature_node) ? feature_node->f.e : NULL;
}

/*
* Return pointer to a face stucture, given the face's name and the
* polyhedron where it lives.
*/
face *findFace(polyhedron *p, char *name)
{
	featureNode *feature_node;

	for (feature_node = p->faces; feature_node && strcmp(name, feature_node->f.f->name); feature_node = feature_node->next);
	return (feature_node) ? feature_node->f.f : NULL;
}

/*
* Return total number of features in a polyhedron
*/
int numFeatures(polyhedron *p)
{
	featureNode *feature_node;
	int n = 0;

	for (feature_node = p->verts; feature_node; feature_node = feature_node->next) n++;
	for (feature_node = p->edges; feature_node; feature_node = feature_node->next) n++;
	for (feature_node = p->faces; feature_node; feature_node = feature_node->next) n++;

	return n;
}

/*
* Return the nth feature of a polyhedron.  The features are numbered starting
* from 0.  Vertices are first, followed by edges, and then faces.
*/
void *nthFeature(polyhedron *p, int n)
{
	featureNode *feature_node;

	for (feature_node = p->verts; n > 0 && feature_node; feature_node = feature_node->next) n--;
	if (feature_node) return feature_node->f.any;
	for (feature_node = p->edges; n > 0 && feature_node; feature_node = feature_node->next) n--;
	if (feature_node) return feature_node->f.any;
	for (feature_node = p->faces; n > 0 && feature_node; feature_node = feature_node->next) n--;
	return feature_node->f.any;
}

/*
* Return a random feature from polyhedron.  There are equal chances that
* a vertex, edge, or face will be returned.  Within each category, all
* features have an equal chance of being selected.
*/
void *randFeat(polyhedron *p)
{
	int type;
	featureNode *fn, *fn2;
	int i, n;

	type = randomNumber(3);
	if (type == 0) fn = p->verts;
	else if (type == 1) fn = p->edges;
	else fn = p->faces;

	// Count number of features in the selected list */
	for (fn2 = fn, n = 0; fn2; n++, fn2 = fn2->next);

	// pick a random feature
	n = randomNumber(n);
	for (i = 0; i < n; i++) fn = fn->next;
	return fn->f.any;
}

/*
* Return radius of a library polyhedron.  The radius of a polyhedron is
* the maximum "norm" of a vertex, where the norm of a vertex is the norm
* of the position vector from the vertex which is rooted at the origin.
*/
float polyhedronRadius(char *name)
{
	float norm, radius;
	featureNode *feature_node;
	int i;

	for (i = 0; i < polyhedronLibraryCount && strcmp(polyhedronLibrary[i].name, name); i++);
	{
		if (i == polyhedronLibraryCount) {
			printf("***** error:  can't find %s in library\n", name);
			exit(1);
		}
	}

	radius = 0.0;
	for (feature_node = polyhedronLibrary[i].verts; feature_node; feature_node = feature_node->next)
	{
		norm = length(feature_node->f.v->coords);
		if (norm > radius) radius = norm;
	}

	return radius;
}

/*
Create a new edge between vertices v1 and v2.  Initialize the edge's
vertex pointers, and it's unit direction vector and length fields.
*/
edge *newEdge(vertex *v1, vertex *v2)
{
	edge *e;

	e = allocEdge;
	e->tag = E;

	if (strcmp(v1->name, v2->name) > 0) SWAP(v1, v2, vertex *);
	e->v1 = v1;
	e->v2 = v2;

	addFeature(e, &v1->edges);
	addFeature(e, &v2->edges);

	e->fl = e->fr = NULL;
	e->cone = NULL;
	e->u = v2->coords - v1->coords;
	e->len = length(e->u);
	e->u = e->u / e->len;
	return e;
}

/*
* Read the polyhedron library from the file whose name is passed in.  The
* library is stored into the global variable polyhedronLibrary.  The global
* variable polyhedronLibraryCount is also initialized.
*/

/*int loadPolyhedronLibrary(char *fname)
{
	FILE *fp;
	char s[80];
	polyhedron *p;
	vertex *v;
	face *f;
	edge *e, *e1, *e2;
	vertex *start, *last;
	featureNode *fn;
	int cont;
	int nameCounter;
	int n;


	printf("\nloading polyhedron library:\n");
	n = 0;
	fopen_s(&fp, fname, "r");

	while (1) {

		do fscanf_s(fp, "%s", s, sizeof(s)); while (!feof(fp) && strcmp(s, "polyhedron"));
		if (feof(fp)) break;

		p = &polyhedronLibrary[n++];
		fscanf_s(fp, "%s", p->name, 20);
		printf_s("%s\n", p->name);
		//mat4copy(mat4I, p->pose);
		p->pose = identity_mat4();
		p->verts = p->edges = p->faces = NULL;

		// read vertices 
		while (1) {
			fscanf_s(fp, "%s", s, sizeof(s));
			if (s[0] == '*') break;
			v = allocVertex;
			v->tag = V;
			strcpy_s(v->name, 20, s);
			v->cone = NULL;
			if (sizeof v->coords.v[0] == sizeof(double))
				fscanf_s(fp, "%f %f %f", &v->coords.v[0], &v->coords.v[1], &v->coords.v[2]);
			else
				fscanf_s(fp, "%f %f %f", &v->coords.v[0], &v->coords.v[1], &v->coords.v[2]);
			v->edges = NULL;
			addFeature(v, &p->verts);
		}

		// read faces 
		nameCounter = 0;
		while (1) {
			fscanf_s(fp, "%s", s, sizeof(s));
			if (s[0] == '*') break;
			if (s[0] == '-') sprintf_s(s, "#F_%d", ++nameCounter);
			f = allocFace;
			f->tag = F;
			strcpy_s(f->name, 20, s);
			f->verts = f->edges = NULL;
			f->cone = NULL;
			addFeature(f, &p->faces);
			getWord(fp, s);
			start = last = findVertex(p, s);
			do {
				if (cont = getWord(fp, s)) v = findVertex(p, s);
				else v = start;
				addFeature(v, &f->verts);
				if (e = findEdge(p, last->name, v->name)) {
					if (e->fl) e->fr = f;
					else e->fl = f;
				}
				else {
					e = newEdge(last, v);
					if (strcmp(last->name, v->name) < 0) e->fl = f;
					else e->fr = f;
					addFeature(e, &p->edges);
				}
				addFeature(e, &f->edges);
				last = v;
			} while (cont);
			// compute face plane coeffs (& outward normal)  
			e1 = f->edges->f.e;
			e2 = f->edges->next->f.e;
			// the "not" below is necessary because we haven't reversed the
			// vertex list to its proper CCW order yet; right now it's in CW order 

			vec3 f_plane;
			if (!(e1->v1 == e2->v2 || e1->v2 == e2->v1))
				f_plane = cross(e1->u, e2->u);
			else f_plane = cross(e2->u, e1->u);
			f_plane = normalise(f_plane);
			f->plane = vec4(f_plane.v[0], f_plane.v[1], f_plane.v[2], f->plane.v[3]);
			f->plane.v[3] = -1 * dot(f->plane, f->verts->f.v->coords);
		}

		// clean up lists 
		for (fn = p->verts; fn; fn = fn->next) reverseFlist(&fn->f.v->edges);
		for (fn = p->faces; fn; fn = fn->next) {
			reverseFlist(&fn->f.f->verts->next);
			reverseFlist(&fn->f.f->edges);
		}
		reverseFlist(&p->verts);
		reverseFlist(&p->edges);
		reverseFlist(&p->faces);

		// build the voronoi regions for the polyhedron 
		buildCones(p);

	}

	fclose(fp);
	polyhedronLibraryCount = n;
	printf("%d polyhedra in library\n\n", n);
	return n;

}*/


/*
* Create polyhedron
*
* Once the polyhedron library has been initialized, we can create instances
* of the various polyhedron in the library.  Any number of instances of the
* same polyhedron may exist, at different positions (poses) in space.  We
* create a polyhedron instance by specifying the name of the reference
* polyhedron in the library (e.g. "cube"), and the name for this particular
* instance (e.g. "cube-1").  The former name is required; the latter one
* is optional, and can be NULL.  A pointer to the newly instantiated polyhedron
* is returned.  The pose matrix is initialized to the identity.  Note that
* all instances of the same library polyhedron share the same vertex, edge,
* and face lists and voronoi structure.
*/
/*polyhedron *createPolyhedron(char *libName, char *name)
{
	polyhedron *newP;
	int i;

	for (i = 0; i < polyhedronLibraryCount &&
		strcmp(polyhedronLibrary[i].name, libName); i++);

	if (i == polyhedronLibraryCount)
	{
		printf("***** error:  can't find %s in library\n", libName);
		return NULL;
	}

	newP = allocPolyhedron;
	strcpy_s(newP->name, 20, name);

	//mat4copy(mat4I, newP->pose);
	newP->pose = identity_mat4();
	newP->verts = polyhedronLibrary[i].verts;
	newP->edges = polyhedronLibrary[i].edges;
	newP->faces = polyhedronLibrary[i].faces;

	return newP;
}*/

polyhedron createPolyhedron(char *name, vector<vec4> bodyVertices, mat4 transformation)
{
	polyhedron *new_polyhedron;
	vertex *v;
	face *f;
	edge *e, *e1, *e2;
	vertex *start, *last;
	featureNode *fn;
	int vertex_count = 0, face_count = 0;

	new_polyhedron = allocPolyhedron;
	strcpy_s(new_polyhedron->name, 20, name);

	new_polyhedron->pose = transformation;
	new_polyhedron->verts = new_polyhedron->edges = new_polyhedron->faces = NULL;

	for (unsigned int i = 0; i < bodyVertices.size(); i++)
	{
		ostringstream oss;
		oss << "v" << vertex_count;
		string vertex_name = oss.str();

		if (!(v = findVertex(new_polyhedron, bodyVertices[i])))
		{
			v = allocVertex;
			v->tag = V;
			strcpy_s(v->name, vertex_name.c_str());
			v->cone = NULL;
			v->coords = bodyVertices[i];
			v->edges = NULL;
			addFeature(v, &new_polyhedron->verts);
			vertex_count++;
		}
	}

	for (unsigned int i = 0; i < (bodyVertices.size() / 3); i++)
	{
		ostringstream oss;
		oss << "f" << face_count;
		string face_name = oss.str();

		f = allocFace;
		f->tag = F;
		strcpy_s(f->name, face_name.c_str());
		f->verts = f->edges = NULL;
		f->cone = NULL;
		addFeature(f, &new_polyhedron->faces);
		face_count++;

		ostringstream vertex_oss;
		vertex_oss << "v" << ((3 * i));
		string vertex_name = vertex_oss.str();

		last = start = findVertex(new_polyhedron, bodyVertices[(3 * i)]);
		//last = start = findVertex(new_polyhedron, _strdup(vertex_name.c_str()));

		for (int j = 1; j <= 3; j++)
		{
			if (j < 3)
			{
				//vertex_oss << "v" << ((3 * i) + j);
				//vertex_name = vertex_oss.str();
				//v = findVertex(new_polyhedron, _strdup(vertex_name.c_str()));
				v = findVertex(new_polyhedron, bodyVertices[(3 * i) + j]);
			}
			else
			{
				v = start;
			}
			addFeature(v, &f->verts);

			if (e = findEdge(new_polyhedron, last->name, v->name))
			{
				if (e->fl) e->fr = f;
				else e->fl = f;
			}
			else
			{
				e = newEdge(last, v);
				if (strcmp(last->name, v->name) < 0) e->fl = f;
				else e->fr = f;
				addFeature(e, &new_polyhedron->edges);
			}
			addFeature(e, &f->edges);
			last = v;
		}

		// Compute face plane coeffs (& outward normal)  
		e1 = f->edges->f.e;
		e2 = f->edges->next->f.e;

		// The "not" below is necessary because we haven't reversed the
		// vertex list to its proper CCW order yet; right now it's in CW order 
		vec3 f_plane;
		
		// TODO: Figure out if the next three lines are incorrect now
		if (!(e1->v1 == e2->v2 || e1->v2 == e2->v1))
			f_plane = cross(e1->u, e2->u);
		else f_plane = cross(e2->u, e1->u);
		
		f_plane = normalise(f_plane);
		f->plane = vec4(f_plane.v[0], f_plane.v[1], f_plane.v[2], f->plane.v[3]);
		f->plane.v[3] = -1 * dot(f->plane, f->verts->f.v->coords);
	}
	
	// Clean up lists 
	for (fn = new_polyhedron->verts; fn; fn = fn->next) reverseFlist(&fn->f.v->edges);

	for (fn = new_polyhedron->faces; fn; fn = fn->next) {
		reverseFlist(&fn->f.f->verts->next);
		reverseFlist(&fn->f.f->edges);
	}

	reverseFlist(&new_polyhedron->verts);
	reverseFlist(&new_polyhedron->edges);
	reverseFlist(&new_polyhedron->faces);

	// Build the voronoi regions for the polyhedron 
	buildCones(new_polyhedron);

	return *new_polyhedron;
}

/*
* Print out all the features of a polyhedron
*/
void dumpPolyhedron(polyhedron *p)
{
	featureNode *fn, *fn2;
	char s[80];

	printf("polyhedron %s ===================\n", p->name);

	printf("vertices:\n");
	for (fn = p->verts; fn; fn = fn->next)
	{
		printf("%-10s  (%+6.2f, %+6.2f, %+6.2f)\n", fn->f.v->name,
			fn->f.v->coords.v[0], fn->f.v->coords.v[1], fn->f.v->coords.v[2]);
		printf("  edges:");
		for (fn2 = fn->f.v->edges; fn2; fn2 = fn2->next)
			printf("  %s.%s", fn2->f.e->v1->name, fn2->f.e->v2->name);
		printf("\n");
	}

	printf("edges:\n");
	for (fn = p->edges; fn; fn = fn->next)
	{
		sprintf_s(s, sizeof(s), "%s.%s", fn->f.e->v1->name, fn->f.e->v2->name);
		printf("%-15s   l=%-10s r=%-10s\n",
			s, fn->f.e->fl->name, fn->f.e->fr->name);
		printf("  l = %+7.2f  u = (%+6.3f, %+6.3f, %+6.3f)\n", fn->f.e->len,
			fn->f.e->u.v[0], fn->f.e->u.v[1], fn->f.e->u.v[2]);
	}

	printf("faces:\n");
	for (fn = p->faces; fn; fn = fn->next)
	{
		printf("%-10s  (%+6.2f, %+6.2f, %+6.2f)  %+6.2f\n", fn->f.f->name,
			fn->f.f->plane.v[0], fn->f.f->plane.v[1],
			fn->f.f->plane.v[2], fn->f.f->plane.v[3]);
		printf("  vertices:");
		for (fn2 = fn->f.f->verts; fn2; fn2 = fn2->next)
			printf("  %s", fn2->f.v->name);
		printf("\n  edges:");
		for (fn2 = fn->f.f->edges; fn2; fn2 = fn2->next)
			printf("  %s.%s", fn2->f.e->v1->name, fn2->f.e->v2->name);
		printf("\n");
	}
}

/*
=============================================================================

Voronoi region routines

These routines are used for building the voronoi region structures
of polyhedra, and for testing for membership of points in these regions.

The three membership routines are vertexConeCheck, edgeConeCheck, and
faceConeCheck, for testing membership of points within vertex, edge, and
face voronoi regions respectively.  Actually, the faceConeCheck only checks
against the side planes of the face's voronoi region; the base plane
must be checked explicitly.  "Points" are simply arrays of 3 floats.

All three membership routines are passed the feature whose cone is being
checked as well as the point.  These routines return true if the point was
in the cone, and false otherwise.  In the latter case, the feature pointer
passed in will be updated with a new feature (of a different type) which
neighboured the old feature, corresponding to the plane which was violated.
We do not update the feature based on the first plane which is found to
be violated, but rather on the plane which is "most violated."  The
amount of violation is measured by the distance of the point from the
plane, assuming the point is on the negative side of the plane.  A point
on the positive side does not violated the plane.  Hence, even if we
find a plane which is violated, we continue to cycle through all planes
of the cone, looking for the one which is most violated.

The voronoi regions are "sticky" in that we make comparisons against
-EPSILON rather than zero.  This provides hysteresis to avoid cycling.
The voronoi planes are oriented so that the "signed distance" of a point
to a voronoi plane is positive if the point is on the region side of the
plane, and negative if the point lies outside the region.

N.B.  The order in which the edge cone planes are built is important!
We add the planes corresponding to vertices (the edge's endpoints) first,
followed by the planes corresponding to the neighboring faces.  This
puts the face planes in the front of the edge's cone list.  We depend
on this order in the edge_face routine when we're checking that the
faces normal lies between the edge's neigbouring faces' normals.

=============================================================================
*/

/*
* Add a plane node to a cone list.  A cone list is a list of planes describing
* the voronoi region of a vertex or an edge, or the side planes of the
* voronoi region of a face.
*/
void addPlane(planeNode *pn, planeNode **cone)
{
	pn->next = *cone;
	*cone = pn;
}

/*
* Flip plane
*
* Planes have an orientation, since we perform checks to see if points lie
* on the proper side of the plane.  flipPlane negates the four parameters
* describing the source plane, to produce a destination plane which is
* identical to the source, except with opposite orientation.
*/
void flipPlane(vec4 src, vec4 &dest)
{
	dest.v[0] = -src.v[0];
	dest.v[1] = -src.v[1];
	dest.v[2] = -src.v[2];
	dest.v[3] = -src.v[3];
}

/*
* Tweak plane normal
*
* Rotate a plane normal epsilon radians about a given axis u.  This
* should be done before the constant term of the plane vector is
* computed, since this term is based on the plane normal.
* u must be a unit vector!
*/
void tweakPlaneNormal(vec3 u, float epsilon, vec3 nOrig, vec4 &nTweak)
{
	mat3 R;
	float cosAngle, sinAngle, v;

	/* compute rotation matrix for rotating 3-vectors an angle eps about u */
	sinAngle = sin(epsilon);
	cosAngle = cos(epsilon);
	v = 1.0f - cosAngle;
	R.m[0] = u.v[0] * u.v[0] * v + cosAngle;
	R.m[4] = u.v[1] * u.v[1] * v + cosAngle;
	R.m[8] = u.v[2] * u.v[2] * v + cosAngle;
	R.m[1] = u.v[0] * u.v[1] * v - u.v[2] * sinAngle;
	R.m[2] = u.v[0] * u.v[2] * v + u.v[1] * sinAngle;
	R.m[3] = u.v[1] * u.v[0] * v + u.v[2] * sinAngle;
	R.m[5] = u.v[1] * u.v[2] * v - u.v[0] * sinAngle;
	R.m[6] = u.v[2] * u.v[0] * v - u.v[1] * sinAngle;
	R.m[7] = u.v[2] * u.v[1] * v + u.v[0] * sinAngle;

	vec3 result = mat3vec3(R, nOrig);
	nTweak.v[0] = result.v[0];
	nTweak.v[1] = result.v[1];
	nTweak.v[2] = result.v[2];
}

/*
* Compute the voronoi region of a vertex.  This creates a list of plane
* nodes (one for each incident edge) and points the vertex's cone field
* to this list.
*/
void computeVertexCone(vertex *v)
{
	planeNode *pn;
	featureNode *fn;
	edge *e;
	vec3 u, tmpV;

	for (fn = v->edges; fn; fn = fn->next)
	{
		e = fn->f.e;
		// Compute vector about which to rotate vertex planes REWRITE
		// this vector points toward the right of edge
		tmpV = e->fl->plane + e->fr->plane;
		u = cross(e->u, tmpV);
		u = normalise(u); // bug fix:  Moritz Breipohl, 21.10.96

		// Construct vertex plane 
		pn = allocPnode;

		if (e->v1 == v)
		{
			tweakPlaneNormal(u, VERTEX_FLARE, e->u, pn->plane);

			pn->plane = pn->plane * -1;
		}
		else
			tweakPlaneNormal(u, -VERTEX_FLARE, e->u, pn->plane);

		pn->plane.v[3] = -1 * dot(v->coords, pn->plane);
		pn->nbr = e;
		addPlane(pn, &v->cone);

		// Construct edge plane
		pn = allocPnode;

		if (e->v1 == v)
			tweakPlaneNormal(u, EDGE_VERTEX_FLARE, e->u, pn->plane);
		else
		{
			tweakPlaneNormal(u, -EDGE_VERTEX_FLARE, e->u, pn->plane);
			pn->plane = pn->plane * -1;
		}

		pn->plane.v[3] = -1 * dot(v->coords, pn->plane);
		pn->nbr = v;
		addPlane(pn, &e->cone);
	}
}

/*
* Compute the voronoi region of a face.  This creates a list of plane
* nodes (one for each edge of the face) and points the faces's cone field
* to this list.  The final plane of the voronoi region is the plane of
* the face itself.  This plane is stored explicitly when the face is
* created.
*/
void computeFaceCone(face *f)
{
	planeNode *pn;
	edge *e;
	vertex *v;
	featureNode *fnE, *fnV;
	vec3 norm;

	// We only do side planes;
	// base plane is already stored explicitly in f->plane 
	for (fnE = f->edges, fnV = f->verts; fnE; fnE = fnE->next, fnV = fnV->next)
	{
		e = fnE->f.e;
		v = fnV->f.v;
		// (e->v1 = v) <==> e goes CCW around f; f is left face of e
		// (e->v2 = v) <==> e goes CW around f; f is right face of e

		// Construct face plane 
		norm = cross(f->plane, e->u);

		pn = allocPnode;

		if (e->v1 == v)
			tweakPlaneNormal(e->u, FACE_FLARE, norm, pn->plane);
		else
		{
			tweakPlaneNormal(e->u, -FACE_FLARE, norm, pn->plane);
			pn->plane = pn->plane * -1;
		}
		pn->plane.v[3] = -1 * dot(v->coords, pn->plane);
		pn->nbr = e;
		addPlane(pn, &f->cone);

		// Construct edge plane 
		pn = allocPnode;
		if (e->v1 == v)
		{
			tweakPlaneNormal(e->u, EDGE_FACE_FLARE, norm, pn->plane);
			pn->plane = pn->plane * -1;
		}
		else
			tweakPlaneNormal(e->u, -EDGE_FACE_FLARE, norm, pn->plane);

		pn->plane.v[3] = -1 * dot(v->coords, pn->plane);
		pn->nbr = f;
		addPlane(pn, &e->cone);
	}
}

/*
* Build the voronoi region structure for an entire polyhedron
*/
void buildCones(polyhedron *p)
{
	featureNode *fn;

	for (fn = p->verts; fn; fn = fn->next)
		computeVertexCone(fn->f.v);
	for (fn = p->faces; fn; fn = fn->next)
		computeFaceCone(fn->f.f);
}

/*
* Print out the voronoi region structure for an entire polyhedron
*/
void dumpCones(polyhedron *p)
{
	featureNode *fn;
	planeNode *pn;

	printf("cone structure of polyhedron %s ===============\n", p->name);

	for (fn = p->verts; fn; fn = fn->next)
	{
		printf("*** vertex %s\n", fn->f.v->name);
		for (pn = fn->f.v->cone; pn; pn = pn->next)
			printf("%+7.3f %+7.3f %+7.3f %+7.3f => %s.%s\n",
				pn->plane.v[0], pn->plane.v[1], pn->plane.v[2], pn->plane.v[3],
				((edge *)pn->nbr)->v1->name,
				((edge *)pn->nbr)->v2->name);
	}

	for (fn = p->edges; fn; fn = fn->next)
	{
		printf("*** edge %s.%s\n", fn->f.e->v1->name, fn->f.e->v2->name);
		for (pn = fn->f.e->cone; pn; pn = pn->next)
			printf("%+7.3f %+7.3f %+7.3f %+7.3f => %s\n",
				pn->plane.v[0], pn->plane.v[1], pn->plane.v[2], pn->plane.v[3],
				(featureTag(pn->nbr) == V) ?
				((VERTEX *)pn->nbr)->name : ((FACE *)pn->nbr)->name);
	}

	for (fn = p->faces; fn; fn = fn->next)
	{
		printf("*** face %s\n", fn->f.f->name);
		for (pn = fn->f.f->cone; pn; pn = pn->next)
		{
			printf("%+7.3f %+7.3f %+7.3f %+7.3f => ",
				pn->plane.v[0], pn->plane.v[1], pn->plane.v[2], pn->plane.v[3]);
			if (pn->nbr) printf("%s.%s\n",
				((edge *)pn->nbr)->v1->name,
				((edge *)pn->nbr)->v2->name);
			else printf("----\n");
		}
	}
}

/*
* Return true if a point lies within vertex v's voronoi region, else false.
* if update is true, update v to a neighboring edge in the latter case.
*/
int vertexConeCheck(vertex **v, vec3 point, int update)
{
	planeNode *pn;
	float min, dot;
	void *feature = NULL;

	min = INFINITY;

	for (pn = (*v)->cone; pn; pn = pn->next)
	{
		if ((dot = planeDist(pn->plane, point)) < min)
		{
			min = dot;
			feature = pn->nbr;
		}
	}

	if (min > -CONE_EPS) return 1;
	// TODO: Check if this is correct
	if (update) *v = (vertex *)feature;
	return 0;
}

/*
* Return true if a point lies within edge e's voronoi region, else false.
* if update is true, update e to a neighboring vertex or face in the latter
* case.
*/
int edgeConeCheck(edge **e, vec3 point, int update)
{
	planeNode *pn;
	float min, dot;
	void *feature = NULL;

	min = INFINITY;
	for (pn = (*e)->cone; pn; pn = pn->next)
	{
		if ((dot = planeDist(pn->plane, point)) < min)
		{
			min = dot;
			feature = pn->nbr;
		}
	}
	if (min >= -CONE_EPS) return 1;
	// TODO: Check if this is correct
	if (update) *e = (edge *)feature;
	return 0;
}

/*
* Return true if a point lies within face f's voronoi region, else false.
* If update is true, update f to a neighboring edge in the latter case.
*/
int faceConeCheck(face **f, vec3 point, int update)
{
	planeNode *pn;
	float min, dot;
	void *feature = NULL;

	min = INFINITY;
	for (pn = (*f)->cone; pn; pn = pn->next)
	{
		if ((dot = planeDist(pn->plane, point)) < min)
		{
			min = dot;
			feature = pn->nbr;
		}
	}

	if (min >= -CONE_EPS) return 1;
	if (update) *f = (face *)feature;
	return 0;
}

/*
=============================================================================

Quick and dirty distance functions and closestToX functions

Occasionally it is necessary to find the closest feature on polyhedron A
to a particular feature on polyhedron B.  We sometimes have to enumerate
over all features on A.  This happens, for instance, when a vertex
on B lies in the "negative voronoi region" of a face on A.

To perform such operations, we use quick and dirty distance functions.
These are seven functions which return the distance between two pairs
of features of a given type.  They are quick and dirty because they
sometimes return infinity rather than an actual distance.  For instance,
the distance between two faces is not even calculated unless the faces
are parallel, for if they are not, we know some some other feature
pair (possibly involving one of the faces) is closer.  Thus Dff returns
infinity if the faces aren't parallel.

These functions are named "Dxy" where x and y are either "v", "e", or
"f", for vertex, edge or face.  Thus, Dev returns the distance between
an edge and a vertex. We use these routines when searching over polyhedron
A for the closest feature to a particular feature on B.  For this reason,
it is most efficient to transform the feature on B to polyhedron A's
frame, rather than transforming every feature on A to B's frame.  In
most of these routines, it is thus assumed that the second feature passed
in (from polyhedron B) has already been transformed to the frame of the
first feature's polyhedron (polyhedron A).  An exception is Dff (see below).

The quick and dirty distance functions are called by the three closestToX
functions:  closestToVert, closestToEdge, and closestToFace.  ClosestToVert
finds the closest feature on A to a vertex on B by first transforming
the vertex on B into A's frame, and making repeated calls to Dvv, Dev, and
Dfv.  Likewise, closestToEdge finds the closest feature on A to an edge
on B by transforming the edge to A's frame and calling Dve, Dee, Dfe while
enumerating over all features of A.  closestToFace works a little
differently, since it is more difficult to transform an entire face into
another frame.  It still enumerates over all features on A, looking for
the closest feature to a face on B, however if the current feature being
examined on A is a vertex or edge, this feature is transformed into B's
frame (rather than transforming the face on B into A's frame).  Hence,
we still make calls to Dfv and Dfe (Dvf and Def don't exist).  If the
current feature on A is a face, we call Dff, which assumes each face is
in its own frame.  This slight inefficiency in the closestToFace case is
no big deal, due to the rarity of calling this function.

Also note in the closestToX functions, we first enumerate over faces, then
edges, then vertices.  This gives priority to features of higher dimension,
which is what we want.  For example, if a face and an edge on A are both
closest to the fixed feature on B, we'll return the face.

=============================================================================
*/

/*
* Return distance between two vertices; the latter has been transformed
* into the frame of the former.
*/
float Dvv(vertex *v1, vertex *v2)
{
	vec3 offset;
	offset = v1->coords - v2->coords;
	return length(offset);
}

/*
* Return distance between an edge and a vertex transformed into the edge's frame.
*/
float Dev(edge *e, vertex *v)
{
	vec3 w, offset;
	float lambda;

	w = v->xcoords - e->v1->coords;
	lambda = dot(w, e->u);
	if (lambda < 0.0)
		lambda = 0.0;
	else if (lambda > e->len)
		lambda = e->len;

	w = displacePoint(e->v1->coords, e->u, lambda);
	offset = v->xcoords - w;

	return length(offset);
}


/*
* Return distance between a face and a vertex transformed into the face's frame,
* or possibly infinity.
*/
float Dfv(face *f, vertex *v)
{
	float dist;
	planeNode *pn;

	for (pn = f->cone; pn; pn = pn->next)
	{
		if (planeDist(pn->plane, v->xcoords) < 0.0) return INFINITY;
	}

	dist = planeDist(f->plane, v->xcoords);
	if (dist < 0.0) return INFINITY;
	return dist;
}


/*
* Return distance between a vertex and an edge transformed into the vertex's frame.
*/
float Dve(vertex *v, edge *e)
{
	vec3 w, offset;
	float lambda;

	w = v->coords - e->v1->xcoords;
	lambda = dot(w, e->xu);
	if (lambda < 0.0)
		lambda = 0.0;
	else if (lambda > e->len)
		lambda = e->len;
	w = displacePoint(e->v1->xcoords, e->xu, lambda);
	offset = v->coords - w;

	return length(offset);
}


/*
* Return distance between two edges; the latter has been transformed into
* the frame of the former.
*/

float Dee(edge *e1, edge *e2)
{
	vec3 cp1, cp2, offset;

	edgeCPs(&e2, &e1, cp2, cp1);
	offset = cp2 - cp1;

	return length(offset);
}


/*
* Return distance between a face and an edge transformed into the face's frame,
* or possibly infinity.
*/
float Dfe(face *f, edge *e)
{
	float h1, h2;
	float min, max;
	vec3 minCut, maxCut;

	if (!polygonCut(f, e, &min, &max)) return INFINITY;
	minCut = displacePoint(e->v1->xcoords, e->xu, min);
	maxCut = displacePoint(e->v1->xcoords, e->xu, max);

	h1 = planeDist(f->plane, minCut);
	h2 = planeDist(f->plane, maxCut);

	if (h1 < 0.0 || h2 < 0.0) return INFINITY;
	return (h1 < h2) ? h1 : h2;
}


/*
* Return distance between two faces, or possibly infinity.  The faces are
* each assumed to be in their own frames, and T12 and T21 are the
* transformation matrices from f1's frame to f2's frame and vice-versa, resp.

* If the faces are parallel, then we return the distance of a vertex on
* one face to the plane of the other.  However, we do this computation both
* ways and take the minimum distance.  For the vertex of a huge face might
* appear to be quite far off the plane of a tiny one due to limited precision.
*/
float Dff(face *f1, face *f2, mat4 T12, mat4 T21)
{
	float k;
	float minDist, dist;
	vec3 w;
	featureNode *fn;
	edge *e;

	w = T21 * f2->plane;
	k = dot(f1->plane, w);
	if (fabs(k) < 1.0 - EPSILON) return INFINITY;

	// At this point we know faces are ||

	minDist = INFINITY;

	// Find closest edge on f1 to f2
	for (fn = f1->edges; fn; fn = fn->next) {
		e = fn->f.e;
		xformEdge(T12, e);
		dist = Dfe(f2, e);
		if (dist < minDist) minDist = dist;
	}

	// now test edges of f2 against f1
	for (fn = f2->edges; fn; fn = fn->next) {
		e = fn->f.e;
		xformEdge(T21, e);
		dist = Dfe(f1, e);
		if (dist < minDist)  minDist = dist;
	}

	return minDist;
}

/*
* Return closest feature on polyhedron p to the vertex v.  Tvp is the
* transformation matrix from v's frame to p's frame.
*/
void *closestToVertex(vertex *v, polyhedron *p, mat4 Tvp)
{
	float dist, minDist;
	void *feature = NULL;
	featureNode *fn;

	// xform vertex to frame of polyhedron p
	xformVertex(Tvp, v);
	minDist = INFINITY;

	for (fn = p->faces; fn; fn = fn->next)
	{
		dist = Dfv(fn->f.f, v);
		if (dist < minDist)
		{
			minDist = dist;
			feature = fn->f.f;
		}
	}
	for (fn = p->edges; fn; fn = fn->next)
	{
		dist = Dev(fn->f.e, v);
		if (dist < minDist)
		{
			minDist = dist;
			feature = fn->f.e;
		}
	}
	for (fn = p->verts; fn; fn = fn->next)
	{
		dist = Dvv(fn->f.v, v);
		if (dist < minDist)
		{
			minDist = dist;
			feature = fn->f.v;
		}
	}
	return feature;
}


/*
* Return closest feature on polyhedron p to the edge e.  Tep is the
* transformation matrix from e's frame to p's frame.
*/
void *closestToEdge(edge *e, polyhedron *p, mat4 Tep)
{
	float dist, minDist;
	void *feature = NULL;
	FNODE *fn;

	// transform edge to frame of polyhedron p
	xformEdge(Tep, e);
	minDist = INFINITY;

	for (fn = p->faces; fn; fn = fn->next)
	{
		dist = Dfe(fn->f.f, e);
		if (dist < minDist)
		{
			minDist = dist;
			feature = fn->f.f;
		}
	}
	for (fn = p->edges; fn; fn = fn->next)
	{
		dist = Dee(fn->f.e, e);
		if (dist < minDist)
		{
			minDist = dist;
			feature = fn->f.e;
		}
	}
	for (fn = p->verts; fn; fn = fn->next)
	{
		dist = Dve(fn->f.v, e);
		if (dist < minDist)
		{
			minDist = dist;
			feature = fn->f.v;
		}
	}
	return feature;
}


/*
* Return closest feature on polyhedron p to the face f.  Tfp is the
* xformation matrix from f's frame to p's frame, and Tpf is the
* inverse transformation.
*/
void *closestToFace(face *f, polyhedron *p, mat4 Tfp, mat4 Tpf)
{
	float dist, minDist;
	void *feature = NULL;
	featureNode *fn;

	minDist = INFINITY;

	for (fn = p->faces; fn; fn = fn->next)
	{
		dist = Dff(fn->f.f, f, Tpf, Tfp);
		if (dist < minDist)
		{
			minDist = dist;
			feature = fn->f.f;
		}
	}
	for (fn = p->edges; fn; fn = fn->next)
	{
		xformEdge(Tpf, fn->f.e);
		dist = Dfe(f, fn->f.e);
		if (dist < minDist)
		{
			minDist = dist;
			feature = fn->f.e;
		}
	}
	for (fn = p->verts; fn; fn = fn->next)
	{
		xformVertex(Tpf, fn->f.v);
		dist = Dfv(f, fn->f.v);
		if (dist < minDist)
		{
			minDist = dist;
			feature = fn->f.v;
		}
	}
	return feature;
}

/*
=============================================================================

Miscellaneous distance support functions

These functions are called by the six basic feature pair checks (i.e.
vertex_vertex, vertex_edge, etc.)  They perform functions such as
computing the closest edge on a face to a given edge, or determining
if two faces overlap.

=============================================================================
*/

/*
* Compute edge closest points
*
* This routine computes closest points between edges e1 and e2, returning
* them in cp1 and cp2.  All possible degenerecies are handled (I think).
* In cases where the closest points are not unique (this can sometimes
* happen with parallel edges), individual closest points are still
* selected.  In these cases, we try to choose sensibly, for instance
* choosing the midpoint of the portion of the edge which is closest to
* the other edge.  We assume e1 has been transformed into e2's frame,
* and both cp1 and cp2 are computed relative to e2's frame.
*/
void edgeCPs(edge **e1, edge **e2, vec3 &cp1, vec3 &cp2)
{
	float k, lambda, lambda1, lambda2;
	vec3 offset, w, y;
	float dot12, dot21;
	float h, t;
	float dist, minDist;

	k = dot((*e1)->xu, (*e2)->u);
	offset = (*e2)->v1->coords - (*e1)->v1->xcoords;

	if (fabs(k) > 1.0 - EPSILON)
	{
		// Case I:  lines || or anti-|| 
		k = (k > 0.0) ? 1.0f : -1.0f;
		lambda1 = lambda2 = INFINITY;
		w = (*e1)->v1->xcoords - (*e2)->v1->coords;
		dot12 = dot(w, (*e2)->u);
		dot21 = -1 * dot(w, (*e1)->xu);

		t = dot12;
		h = dot12 + k * (*e1)->len;

		if (t <= 0.0)
		{
			if (h <= 0.0) lambda2 = 0.0;
			else if (h >= (*e2)->len) lambda2 = (*e2)->len / 2;
		}
		else if (t >= (*e2)->len)
		{
			if (h <= 0.0) lambda2 = (*e2)->len / 2;
			else if (h >= (*e2)->len) lambda2 = (*e2)->len;
		}

		t = dot21;
		h = dot21 + k * (*e2)->len;

		if (t <= 0.0)
		{
			if (h <= 0.0) lambda1 = 0.0;
			else if (h >= (*e1)->len) lambda1 = (*e1)->len / 2;
		}
		else if (t >= (*e1)->len)
		{
			if (h <= 0.0) lambda1 = (*e1)->len / 2;
			else if (h >= (*e1)->len) lambda1 = (*e1)->len;
		}

		if (lambda1 != (float)INFINITY && lambda2 == (float)INFINITY)
			lambda2 = dot12 + k * lambda1;
		else if (lambda2 != (float)INFINITY && lambda1 == (float)INFINITY)
			lambda1 = dot21 + k * lambda2;
		else if (lambda1 == (float)INFINITY)
		{
			// N.B. lambda2 also == INFINITY
			if (t < 0.0) lambda1 = h / 2;
			else if (t >(*e1)->len) lambda1 = (h + (*e1)->len) / 2;
			else if (h < 0.0) lambda1 = t / 2;
			else lambda1 = (t + (*e1)->len) / 2;
			lambda2 = dot12 + k * lambda1;
		}

		// Compute cps based on lambdas

		if (lambda1 <= 0.0) cp1 = (*e1)->v1->xcoords;
		else if (lambda1 >= (*e1)->len) cp1 = (*e1)->v2->xcoords;
		else cp1 = displacePoint((*e1)->v1->xcoords, (*e1)->xu, lambda1);

		if (lambda2 <= 0.0) cp2 = (*e2)->v1->coords;
		else if (lambda2 >= (*e2)->len) cp2 = (*e2)->v2->coords;
		else cp2 = displacePoint((*e2)->v1->coords, (*e2)->u, lambda2);
	}
	else
	{
		// Case II:  lines not || 
		w = (*e2)->u * k;
		w = (*e1)->xu - w;
		lambda1 = dot(w, offset) / (1 - k * k);
		w = (*e1)->xu * k;
		w = (*e2)->u - w;
		lambda2 = -1 * dot(w, offset) / (1 - k * k);

		if (lambda1 >= 0.0 && lambda1 <= (*e1)->len &&
			lambda2 >= 0.0 && lambda2 <= (*e2)->len)
		{
			cp1 = displacePoint((*e1)->v1->xcoords, (*e1)->xu, lambda1);
			cp2 = displacePoint((*e2)->v1->coords, (*e2)->u, lambda2);
			return;
		}

		// (lambda1, lambda2) not in allowable region; check boundaries 

		minDist = INFINITY;
		if (lambda1 < 0.0)
		{
			// Check boundary:  lambda1 = 0 
			w = (*e1)->v1->xcoords - (*e2)->v1->coords;
			lambda = dot(w, (*e2)->u);
			if (lambda < 0.0) lambda = 0.0;
			else if (lambda >(*e2)->len) lambda = (*e2)->len;
			w = displacePoint((*e2)->v1->coords, (*e2)->u, lambda);
			y = (*e1)->v1->xcoords - w;
			dist = length(y);
			if (dist < minDist)
			{
				minDist = dist;
				cp1 = (*e1)->v1->xcoords;
				cp2 = w;
			}
		}
		else if (lambda1 >(*e1)->len)
		{
			// Check boundary:  lambda1 = length1 
			w = (*e1)->v2->xcoords - (*e2)->v1->coords;
			lambda = dot(w, (*e2)->u);
			if (lambda < 0.0) lambda = 0.0;
			else if (lambda >(*e2)->len) lambda = (*e2)->len;
			w = displacePoint((*e2)->v1->coords, (*e2)->u, lambda);
			y = (*e1)->v2->xcoords - w;
			dist = length(y);
			if (dist < minDist) {
				minDist = dist;
				cp1 = (*e1)->v2->xcoords;
				cp2 = w;
			}
		}

		if (lambda2 < 0.0) {
			/* check boundary:  lambda2 = 0 */
			w = (*e2)->v1->coords - (*e1)->v1->xcoords;
			lambda = dot(w, (*e1)->xu);
			if (lambda < 0.0) lambda = 0.0;
			else if (lambda >(*e1)->len) lambda = (*e1)->len;
			w = displacePoint((*e1)->v1->xcoords, (*e1)->xu, lambda);
			y = (*e2)->v1->coords - w;
			dist = length(y);
			if (dist < minDist) {
				minDist = dist;
				cp2 = (*e2)->v1->coords;
				cp1 = w;
			}
		}
		else if (lambda2 >(*e2)->len) {
			/* check boundary:  lambda2 = length2 */
			w = (*e2)->v2->coords - (*e1)->v1->xcoords;
			lambda = dot(w, (*e1)->xu);
			if (lambda < 0.0) lambda = 0.0;
			else if (lambda >(*e1)->len) lambda = (*e1)->len;
			w = displacePoint((*e1)->v1->xcoords, (*e1)->xu, lambda);
			y = (*e2)->v2->coords - w;
			dist = length(y);
			if (dist < minDist)
			{
				minDist = dist;
				cp2 = (*e2)->v2->coords;
				cp1 = w;
			}
		}
	}
}

/*
* Closest edge or vertex on face
*
* Return the edge or vertex of face f which is closest to edge e.
* We assume e has been transformed to f's frame.
*/
void *closestEdgeOrVertexOnFace(face *f, edge *e)
{
	featureNode *fn;
	float dist, minDist;
	edge *e2, *closest = NULL;
	vec3 cp1, cp2, tmp;

	minDist = INFINITY;
	for (fn = f->edges; fn; fn = fn->next)
	{
		e2 = fn->f.e;
		edgeCPs(&e, &e2, cp1, cp2);
		tmp = cp2 - cp1;
		dist = length(tmp);
		if (dist < minDist)
		{
			minDist = dist;
			if (vec3equal(cp2, e2->v1->coords)) closest = (edge *)e2->v1;
			else if (vec3equal(cp2, e2->v2->coords)) closest = (edge *)e2->v2;
			else closest = e2;
		}
	}
	return closest;
}

/*
* Find the closest pair of edges between faces f1 and f2; return them
* in closest1 & closest2.  T12 is the transformation matrix from f1's
* frame to f2's frame.  This is just a dumb n^2 algorithm which loops
* through all pairs of edges.

* A minor complication is that we favor edge pairs for which the closest
* points aren't endpoints of the edge.  So in addition to computing the
* distance between two edges, we compute the "subdistance," an integer
* 0-2 which indicates how many of the closest points between the edge
* pair are endpoints.  When looking for the closest pair, we first
* check the distance of a pair to the minimum distance found thus far.  If
* these distances are close, we check the subdistances.
*/

void closestEdges(face *f1, face *f2, mat4 T12, edge **closest1, edge **closest2)
{
	featureNode *fn1, *fn2;
	float dist, minDist;
	int subDist, minSubDist;
	edge *e1, *e2;
	vec3 cp1, cp2, tmp;

	minDist = INFINITY;
	for (fn1 = f1->edges; fn1; fn1 = fn1->next)
	{
		e1 = fn1->f.e;
		xformEdge(T12, e1);
		for (fn2 = f2->edges; fn2; fn2 = fn2->next)
		{
			e2 = fn2->f.e;
			edgeCPs(&e1, &e2, cp1, cp2);
			tmp = cp2 - cp1;
			dist = length(tmp);
			// Adjustment to favor edges for which closest points aren't endpoints
			subDist = (vec3equal(cp1, e1->v1->xcoords) || vec3equal(cp1, e1->v2->xcoords)) ? 1 : 0;

			if (vec3equal(cp2, e2->v1->coords) || vec3equal(cp2, e2->v2->coords))
				subDist++;
			if ((dist < minDist) || (dist < minDist + EPSILON && subDist < minSubDist))
			{
				minDist = dist;
				minSubDist = subDist;
				*closest1 = e1;
				*closest2 = e2;
			}
		}
	}
}

/*
* Scan over the edges and vertices of face 2, and return the one which
* is closest to the plane of f1.  T12 & T21 are the usual transformation
* matrices.  Edges are favored over vertices when they are parallel to
* the face.
*/
void *closestToFacePlane(face *f1, face *f2, mat4 T12, mat4 T21)
{
	FNODE *fn;
	edge *e;
	float dist, minDist, dist1, dist2;
	void *feature, *closest = NULL;
	vec3 norm;

	minDist = INFINITY;
	xformVector(T12, f1->plane, norm);
	for (fn = f2->edges; fn; fn = fn->next)
	{
		e = fn->f.e;

		xformVertex(T21, e->v1);
		dist1 = fabs(planeDist(f1->plane, e->v1->xcoords));
		xformVertex(T21, e->v2);
		dist2 = fabs(planeDist(f1->plane, e->v2->xcoords));
		if (dist1 < dist2 - EPSILON)
		{
			// v1 is closer 
			dist = dist1;
			feature = e->v1;
		}
		else if (dist2 < dist1 - EPSILON)
		{
			// v2 is closer 
			dist = dist2;
			feature = e->v2;
		}
		else
		{
			// Edge is parallel to face 
			dist = (dist1 < dist2) ? dist1 : dist2;
			feature = e;
		}

		if ((dist < minDist) || (dist < minDist + EPSILON &&
			featureTag(closest) == V && featureTag(feature) == E))
		{
			minDist = dist;
			closest = feature;
		}
	}
	return closest;
}

/*
* Determine if edge e intersects face f's cone.  Return true if it
* does, else false.  If true, min and max are parameter values which
* bound the segment of e contained in the cone.  Note that
* 0 <= min < max <= length(e).  We assume e has been transformed to
* f's frame.
*/
int polygonCut(face *f, edge *e, float *min, float *max)
{
	planeNode *pn;
	float lambda;
	float denom;

	*min = 0.0;
	*max = e->len;

	for (pn = f->cone; pn; pn = pn->next)
	{
		denom = dot(e->xu, pn->plane);
		if (denom > 0.0)
		{
			lambda = -planeDist(pn->plane, e->v1->xcoords) / denom;
			if (lambda > *min)
			{
				*min = lambda;
				if (*min > *max) return 0;
			}
		}
		else if (denom < 0.0)
		{
			lambda = -planeDist(pn->plane, e->v1->xcoords) / denom;
			if (lambda < *max)
			{
				*max = lambda;
				if (*max < *min) return 0;
			}
		}
		else
		{
			// denom = 0.0 ; lines are parallel
			if (planeDist(pn->plane, e->v1->xcoords) < 0.0) return 0;
		}
	}
	return 1;
}

/*
* Return true if faces overlap, else false.  This is only called with
* pairs of parallel faces, and by "overlap" we mean the projection of
* one face onto the plane of the other overlaps the other face.
*
* We return true iff:
* some edge of f1 intersects the cone of f2, OR
* an abritrarily chosen vertex of f2 lies in the cone of f1
*/
int faceOverlap(face *f1, face *f2, mat4 T12, mat4 T21)
{
	featureNode *fn;
	float min, max;

	for (fn = f1->edges; fn; fn = fn->next)
	{
		xformEdge(T12, fn->f.e);
		if (polygonCut(f2, fn->f.e, &min, &max)) return 1;
	}

	// No dice - maybe f2 lies completely within f1 
	xformVertex(T21, f2->verts->f.v);
	return faceConeCheck(&f1, f2->verts->f.v->xcoords, 0);
}

/*
=============================================================================

feature pair checks and closest feature routine

The six feature pair check routines test if a given pair of features
is the closest pair of features between two polyhedra.  Each one handles
a different combination of feature types.  The six routines are:
vertex_vertex, vertex_edge, vertex_face, edge_edge, edge_face, and face_face.
Each of these routines takes as parameters the two features as well as
two transformation matrices for transforming features from one frame to
the other and vice versa.  T12 is always the transformation matrix from
the frame of the first feature in the parameter list to the frame of
the second feature in the parameter list; T21 is the inverse transformation.
The check routines which involve faces also take pointer(s) to the
polyhedron(a) of the face(s), since it is necessary to scan over the entire
polyhedron when a face's base plane is violated by a feature from the other
polyhedron.

If pair of features passed into a check routine is indeed a closest
feature pair, the distance between these features is returned.  If the
features are not a closest feature pair, -INFINITY is returned, and
one or possibly both of the closest features are updated to form a
new pair which is guaranteed to be closer than the previous pair.

The routine closestFeatures is the main loop of the distance algorithm.
It iteratively "walks" around both polyhedra, calling the appropriate
feature pair check routines, until it finds the closest pair of features.
It then returns this pair.

=============================================================================
*/

/*
* check if two vertices are a closest feature pair
*/
float vertex_vertex(vertex **v1, vertex **v2, mat4 T12, mat4 T21)
{
	vec3 offset;

	// Check if v1 lies in v2's cone
	xformVertex(T12, *v1);
	if (!vertexConeCheck(v2, (*v1)->xcoords, 1)) return -INFINITY;

	// Check if v2 lies in v1's cone 
	xformVertex(T21, *v2);
	if (!vertexConeCheck(v1, (*v2)->xcoords, 1)) return -INFINITY;

	offset = (*v1)->coords - (*v2)->xcoords;
	return length(offset);
}

/*
* Check if a vertex and an edge are a closest feature pair
*/
float vertex_edge(vertex **v, edge **e, mat4 T12, mat4 T21)
{
	vec3 tmp, tmp2, offset;
	float lambda;

	// check if v lies in e's cone
	xformVertex(T12, *v);
	if (!edgeConeCheck(e, (*v)->xcoords, 1)) return -INFINITY;;

	// Compute closest point on e to v
	tmp = (*v)->xcoords - (*e)->v1->coords;
	lambda = dot(tmp, (*e)->u);
	if (lambda < 0.0) lambda = 0.0;
	else if (lambda >(*e)->len) lambda = (*e)->len;
	tmp = displacePoint((*e)->v1->coords, (*e)->u, lambda);

	// Check if closest point on e lies in v's cone
	xformPoint(T21, tmp, tmp2);
	if (!vertexConeCheck(v, tmp2, 1)) return -INFINITY;
	offset = (*v)->coords - tmp2;
	return length(offset);
}

/*
* Check if two edges are a closest feature pair
*/
float edge_edge(edge **e1, edge **e2, mat4 T12, mat4 T21)
{
	vec3 cp1, cp2, offset, xformed;

	xformEdge(T12, *e1);
	edgeCPs(e1, e2, cp1, cp2);

	// Check if closest point on e1 lies in e2's cone
	if (!edgeConeCheck(e2, cp1, 1)) return -INFINITY;

	// Check if closest point on e2 lies in e1's cone 
	xformPoint(T21, cp2, xformed);
	if (!edgeConeCheck(e1, xformed, 1)) return -INFINITY;

	offset = cp1 - cp2;
	return length(offset);
}

/*
* Check if a vertex and a face are a closest feature pair
*/
float vertex_face(vertex **v, face **f, mat4 T12, mat4 T21, polyhedron *facePoly)
{
	vec3 tmp, tmp2;
	float dist;

	// check if v lies in f's cone
	xformVertex(T12, *v);
	if (!faceConeCheck(f, (*v)->xcoords, 1)) return -INFINITY;
	dist = planeDist((*f)->plane, (*v)->xcoords);
	if (dist < 0.0)
	{
		*f = (face *)closestToVertex(*v, facePoly, T12);
		return -INFINITY;
	}

	// 	Compute closest point on f to v, i.e. proj. of v in f's plane.
	// At this point, we know the projected point lies within f
	tmp = (*f)->plane * dist;
	tmp = (*v)->xcoords - tmp;

	// check if closest point on f lies in v's cone
	xformPoint(T21, tmp, tmp2);
	if (!vertexConeCheck(v, tmp2, 1)) return -INFINITY;

	return dist;
}

float edge_face(edge **e, face **f, mat4 T12, mat4 T21, polyhedron *facePoly)
{
	planeNode *pn;
	vec3 w, tmp;
	float h1, h2, absH1, absH2;
	float epsilon;
	vec3 tmpV1, tmpV2;
	float min, max;

	// Compute distance of v1 and v2 from face
	xformEdge(T12, *e);
	h1 = planeDist((*f)->plane, (*e)->v1->xcoords);
	h2 = planeDist((*f)->plane, (*e)->v2->xcoords);
	absH1 = fabs(h1);
	absH2 = fabs(h2);

	// Compute a variable epsilon which will prevent
	// edge-face / vertex-face cycling  and edge-face edge-edge cycling
	// (1.5 is a safety factor)
	if (EDGE_FACE_FLARE > EDGE_VERTEX_FLARE /*VERT_FLARE*/)
		epsilon = /*1.5*/ (float)1.0e-4 + (*e)->len * (float)sin(EDGE_FACE_FLARE);
	else epsilon = /*1.5*/ (float)1.0e-4 + (*e)->len * (float)sin(EDGE_VERTEX_FLARE /*VERT_FLARE*/);

	// Fix for edges that would otherwise look || to f if h1 = -h2
	if (((h1 < 0.0 && h2 > 0.0) || (h1 > 0.0 && h2 < 0.0))
		&& (fabs(h1 - h2) > epsilon))
	{
		*f = (face *)closestEdgeOrVertexOnFace(*f, *e);
		return -INFINITY;
	}

	if (absH1 < absH2 - epsilon)
	{
		// Case I: v1 is closer 
		// Check if v1 lies in f's cone 
		for (pn = (*f)->cone; pn; pn = pn->next)
		{
			if (planeDist(pn->plane, (*e)->v1->xcoords) < 0.0) break;
		}
		if (pn) *f = (face *)closestEdgeOrVertexOnFace(*f, *e);
		else if (h1 < 0.0) *f = (face *)closestToEdge(*e, facePoly, T12);
		else *e = (edge *)(*e)->v1;
		return -INFINITY;
	}

	if (absH2 < absH1 - epsilon)
	{
		// Case II: v2 is closer
		// Check if v2 lies in f's cone
		for (pn = (*f)->cone; pn; pn = pn->next)
		{
			if (planeDist(pn->plane, (*e)->v2->xcoords) < 0.0) break;
		}
		if (pn) *f = (face *)closestEdgeOrVertexOnFace(*f, *e);
		else if (h2 < 0.0) *f = (face *)closestToEdge(*e, facePoly, T12);
		else *e = (edge *)(*e)->v2;
		return -INFINITY;
	}

	else
	{
		// Case III: e is parallel to face 

		if (!polygonCut(*f, *e, &min, &max))
		{
			*f = (face *)closestEdgeOrVertexOnFace(*f, *e);
			return -INFINITY;
		}

		tmpV1 = displacePoint((*e)->v1->xcoords, (*e)->xu, min);
		h1 = planeDist((*f)->plane, tmpV1);
		tmpV2 = displacePoint((*e)->v1->xcoords, (*e)->xu, max);
		h2 = planeDist((*f)->plane, tmpV2);
		if (h1 < 0.0 || h2 < 0.0)
		{
			*f = (face *)closestToEdge(*e, facePoly, T12);
			return -INFINITY;
		}

		// With flared planes, Ming's algo must be changed.  We must compare the
		// face's normal directly to the two face plane normals of the edge's cone,
		// rather than comparing it to a vector derived from the edge direction and
		// the normals of the neighboring faces.
		// Note we assume the face planes are the first two in the edge's cone list!
		xformVector(T21, (*f)->plane, w);
		if (dot(w, (*e)->cone->plane) > 0.0)
		{
			*e = (edge *)(*e)->cone->nbr;
			return -INFINITY;
		}
		if (dot(w, (*e)->cone->next->plane) > 0.0)
		{
			*e = (edge *)(*e)->cone->next->nbr;
			return -INFINITY;
		}

		if (h1 <= h2)
		{
			tmp = (*f)->plane * h1;
			return h1;
		}
		else
		{
			tmp = (*f)->plane * h2;
			return h2;
		}
	}
}

/*
* Check if two faces are a closest feature pair
*
* The parallel, overlapping face-face case is rather complicated.
* We first try to compute
* the distance from face A to face B by finding the smallest positive dist of
* any vertex of A from the plane of B or any vertex of B from the plane of A.
* Negative distances are ignored under the assumption that they are just due
* to numerical errors that result when the vertex of A is very laterally
* distant from B (e.g. inifinite planes), unless all distances are negative,
* in which case we assume faces have interpenetrated.
*/
float face_face(face **f1, face **f2, mat4 T12, mat4 T21, polyhedron *face1poly, polyhedron *face2poly)
{
	planeNode *pn;
	featureNode *fn;
	edge *e1, *e2, *e;
	vec3 w;
	float k;
	void *feature = NULL;
	float min, max;
	int penetration;
	float dist, minDist;
	vec3 tmpV, tmpV2;

	xformVector(T21, (*f2)->plane, w);
	k = dot((*f1)->plane, w);
	if (fabs(k) > 1.0 - EPSILON)
	{
		// Case I: faces are parallel
		if (faceOverlap(*f1, *f2, T12, T21))
		{
			// Modified this case.  Used to determine if plane A was on (parallel)
			// plane B's inner side if the first vertex of A was on B's negative
			// side.  This cause probs when A is very large compared to B - even
			// if planes are roughly parallel, a distant vertex of A could appear
			// to lie quite far into B's negative side.

			// Fix: we clip one polygon to area of other one, therby
			// avoiding checks of vertices which are located far away from
			// one of the polygons. 

			minDist = INFINITY;
			penetration = 0;

			// Look at where f2's cone cuts f1's edges
			for (fn = (*f1)->edges; fn; fn = fn->next)
			{
				e = fn->f.e;
				xformEdge(T12, e);
				if (polygonCut(*f2, e, &min, &max))
				{
					tmpV = displacePoint(e->v1->xcoords, e->xu, min);
					dist = planeDist((*f2)->plane, tmpV);
					if (dist < minDist)
					{
						if (penetration = (dist < 0.0)) break;
						minDist = dist;
						tmpV2 = (*f2)->plane * dist;
					}
					tmpV = displacePoint(e->v1->xcoords, e->xu, max);
					dist = planeDist((*f2)->plane, tmpV);
					if (dist < minDist)
					{
						if (penetration = (dist < 0.0)) break;
						minDist = dist;
						tmpV2 = (*f2)->plane * dist;
					}
				}
			}
			if (penetration)
			{
				*f2 = (face *)closestToFace(*f1, face2poly, T12, T21);
				return -INFINITY;
			}

			// Look at where f1's cone cuts f2's edges
			for (fn = (*f2)->edges; fn; fn = fn->next)
			{
				e = fn->f.e;
				xformEdge(T21, e);
				if (polygonCut(*f1, e, &min, &max))
				{
					tmpV = displacePoint(e->v1->xcoords, e->xu, min);
					dist = planeDist((*f1)->plane, tmpV);
					if (dist < minDist)
					{
						if (penetration = (dist < 0.0)) break;
						minDist = dist;
						tmpV2 = (*f1)->plane * dist;
					}
					tmpV = displacePoint(e->v1->xcoords, e->xu, max);
					dist = planeDist((*f1)->plane, tmpV);
					if (dist <  minDist) {
						if (penetration = (dist < 0.0)) break;
						minDist = dist;
						tmpV2 = (*f1)->plane * dist;
					}
				}
			}
			if (penetration)
			{
				*f1 = (face *)closestToFace(*f2, face1poly, T21, T12);
				return -INFINITY;
			}

			return minDist;
		}

		closestEdges(*f1, *f2, T12, &e1, &e2);
		*f1 = (face *)e1;
		*f2 = (face *)e2;
		return -INFINITY;
	}

	else
	{
		// Case II: faces not parallel

		// Return f2 & edge/vert of f1
		feature = closestToFacePlane(*f2, *f1, T21, T12);
		if (featureTag(feature) == V)
		{
			xformVertex(T12, (vertex *)feature);
			for (pn = (*f2)->cone; pn; pn = pn->next)
			{
				if (planeDist(pn->plane, ((vertex *)feature)->xcoords) < 0.0) break;
			}
			if (!pn)
			{
				*f1 = (face *)feature;
				return -INFINITY;
			}
		}
		else
		{
			e = (EDGE *)feature;
			xformEdge(T12, e);
			/*
			h1 = planeDist((*f2)->plane, e->v1->xcoords);
			h2 = planeDist((*f2)->plane, e->v2->xcoords);
			vectScale((*f2)->plane, tmp, h1);
			vectSub(e->v1->xcoords, tmp, w1);
			vectScale((*f2)->plane, tmp, h2);
			vectSub(e->v2->xcoords, tmp, w2);
			if (polygonCut(*f2, w1, w2)) {
			*/
			if (polygonCut(*f2, e, &min, &max))
			{
				*f1 = (face *)feature;
				return -INFINITY;
			}
		}

		/* return f1 & edge/vert of f2 */
		feature = closestToFacePlane(*f1, *f2, T12, T21);
		if (featureTag(feature) == V)
		{
			xformVertex(T21, (vertex *)feature);
			for (pn = (*f1)->cone; pn; pn = pn->next)
			{
				if (planeDist(pn->plane, ((vertex *)feature)->xcoords) < 0.0) break;
			}
			if (!pn)
			{
				*f2 = (face *)feature;
				return -INFINITY;
			}
		}
		else
		{
			e = (edge *)feature;
			xformEdge(T21, e);
			/*
			h1 = planeDist((*f1)->plane, e->v1->xcoords);
			h2 = planeDist((*f1)->plane, e->v2->xcoords);
			vectScale((*f1)->plane, tmp, h1);
			vectSub(e->v1->xcoords, tmp, w1);
			vectScale((*f1)->plane, tmp, h2);
			vectSub(e->v2->xcoords, tmp, w2);
			if (polygonCut(*f1, w1, w2)) {
			*/
			if (polygonCut(*f1, e, &min, &max)) {
				*f2 = (face *)feature;
				return -INFINITY;
			}
		}
		closestEdges(*f1, *f2, T12, &e1, &e2);
		*f1 = (face *)e1;
		*f2 = (face *)e2;
		return -INFINITY;
	}
}

/*
* Try to get a distance between features, even if the closest
* feature algorithm is cycling.
*
* After the cycle detector is reset (by calling with reset = 1),
* cycleDetector() will repeatedly return -INFINITY, until a cycle
* is detected, at which time it will return a positive distance
* - posisbly infinity - according to the smallest distance found so far.
* A cycle is detected when the distance between the current features
* fails to be lower than the minimum dist found thus far at least 5 times.
*/
float cycleDetector(void **feat1, void **feat2, mat4 T12, mat4 T21, int reset)
{
	float dist;
	int flag;
	char name1[20], name2[20];

	static int cycleCounter;
	static float minDist;
	static void *bestCF1, *bestCF2;

	if (reset)
	{
		minDist = INFINITY;
		cycleCounter = 0;
		return 0.0;  // Return value irrelevant here
	}

	featureName(*feat1, name1);
	featureName(*feat2, name2);
	printf("cycle? : feat1 = %-20s feat2 = %-s\n", name1, name2);

	flag = (featureTag(*feat1) << 2) + featureTag(*feat2);
	switch (flag)
	{
	case (V << 2) + V:
		xformVertex(T21, (vertex *)feat2);
		dist = Dvv((vertex *)feat1, (vertex *)feat2);
		break;
	case (V << 2) + E:
		xformEdge(T21, (edge *)feat2);
		dist = Dve((vertex *)feat1, (edge *)feat2);
		break;
	case (V << 2) + F:
		xformVertex(T12, (vertex *)feat1);
		dist = Dfv((face *)feat2, (vertex *)feat1);
		break;
	case (E << 2) + V:
		xformVertex(T21, (vertex *)feat2);
		dist = Dev((edge *)feat1, (vertex *)feat2);
		break;
	case (E << 2) + E:
		xformEdge(T21, (edge *)feat2);
		dist = Dee((edge *)feat1, (edge *)feat2);
		break;
	case (E << 2) + F:
		xformEdge(T12, (edge *)feat1);
		dist = Dfe((face *)feat2, (edge *)feat1);
		break;
	case (F << 2) + V:
		xformVertex(T21, (vertex *)feat2);
		dist = Dfv((face *)feat1, (vertex *)feat2);
		break;
	case (F << 2) + E:
		xformEdge(T21, (edge *)feat2);
		dist = Dfe((face *)feat1, (edge *)feat2);
		break;
	case (F << 2) + F:
		dist = Dff((face *)feat1, (face *)feat2, T12, T21);
		break;
	}

	printf("dist=%+14.10f min=%+14.10f\n", dist, minDist);
	if (dist >= minDist) cycleCounter++;
	else
	{
		minDist = dist;
		bestCF1 = *feat1;
		bestCF2 = *feat2;
	}

	if (cycleCounter == 5)
	{
		// Cycle detected!
		*feat1 = bestCF1;
		*feat2 = bestCF2;

		// Uncomment next line to find out how often cycles are occuring 
		// printf("cycle detector returns %+12.6f\n", minDist); 

		return minDist;
	}

	else return -INFINITY;
}

/*
* Closest feature pair
*
* Compute the closest feature pair for polyhedra poly1 & poly2.  The initial
* feature pair passed is passed in via feat1 & feat2, and upon termination
* these are updated to reflect the new closest feature pair.  The routine
* returns the distance between the pair of closest features.
*
* We keep track of the number of iterations of the algorithm (i.e. how
* many "steps" we take as we walk around the polyhedra, looking for
* closest features).  If this count reaches CYCLE_CNT, we assume we
* are stuck in a cycle, and invoke the cycle detection system.
* Under normal circumstances, this should rarely occur.  If the polyhedra
* have high numbers of facets (tens of thousands) or experience large
* changes in pose between calls to closestFeatures, CYCLE_CNT may have
* to be increased.
*/
#define CYCLE_CNT 500

float closestFeatures(polyhedron *poly1, void **feat1, polyhedron *poly2, void **feat2)
{
	int type1, type2;
	float dist;
	mat4 T12, T21, inv;
	int cycleChk;
	char name1[60], name2[60];

	// Compute transformation matrices between the two polyhedron frames
	matInvXform(poly2->pose, inv);
	matMultXform(inv, poly1->pose, T12);
	matInvXform(T12, T21);

	cycleChk = 0;

	do {

		if (cycleChk++ == CYCLE_CNT) cycleDetector(feat1, feat2, T12, T21, 1);
		else if (cycleChk > CYCLE_CNT)
		{
			dist = cycleDetector(feat1, feat2, T12, T21, 0);
			if (dist > 0.0) return dist;
		}

		type1 = featureTag(*feat1);
		type2 = featureTag(*feat2);

		featureName(*feat1, name1);
		featureName(*feat2, name2);
		printf("cycle? : feat1 = %s feat2 = %s\n", name1, name2);

		if (type1 == V)
		{
			if (type2 == V) dist = vertex_vertex((vertex **)feat1, (vertex **)feat2, T12, T21);
			else if (type2 == E) dist = vertex_edge((vertex **)feat1, (edge **)feat2, T12, T21);
			else dist = vertex_face((vertex **)feat1, (face **)feat2, T12, T21, poly2);
		}
		else if (type1 == E)
		{
			if (type2 == V) dist = vertex_edge((vertex **)feat2, (edge **)feat1, T21, T12);
			else if (type2 == E) dist = edge_edge((edge **)feat1, (edge **)feat2, T12, T21);
			else dist = edge_face((edge **)feat1, (face **)feat2, T12, T21, poly2);
		}
		else
		{
			if (type2 == V) dist = vertex_face((vertex **)feat2, (face **)feat1, T21, T12, poly1);
			else if (type2 == E)
				dist = edge_face((edge **)feat2, (face **)feat1, T21, T12, poly1);
			else dist = face_face((face **)feat1, (face **)feat2, T12, T21, poly1, poly2);
		}

	} while (dist < 0.0);

	return dist;
}

/*
* Intialize closest feature pair
*
* Functionally, this is identical to closestFeatures, except we call this
* routine the first time, when we have not yet initialized feat1 and feat2.
* This routine simply chooses a starting feature on each of the polyhedra,
* and then calls closestFeatures.
*
* We choose the starting features in a fairly dumb manner, picking the
* first vertex in each polyhedron's vertex list.  This is totally arbitrary.
* There may be smarter ways of doing it, but it really doesn't matter, since
* we only perform this operation once, and then all subsequent calls are to
* closestFeatures.
*/
float closestFeaturesInit(polyhedron *poly1, void **feat1, polyhedron *poly2, void **feat2)
{
	*feat1 = poly1->verts->f.v;
	*feat2 = poly2->verts->f.v;
	return closestFeatures(poly1, feat1, poly2, feat2);
}