#include <string.h>
#include <stdlib.h>

#include "ClosestFeatures.h"

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

	do fscanf(fp, "%c", &c); while ((c == ' ') || c == '\t');
	if (c == '\n') return 0;
	ungetc(c, fp);
	fscanf(fp, "%s", s);
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
	if (type == V) strcpy(name, ((vertex *)feature)->name);
	else if (type == F) strcpy(name, ((face *)feature)->name);
	else sprintf(name, "%s.%s", ((edge *)feature)->v1->name, ((edge *)feature)->v2->name);
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
int numFeats(polyhedron *p)
{
	featureNode *feature_node;
	int n = 0;

	for (feature_node = p->verts; feature_node; feature_node = feature_node->next) n++;
	for (feature_node = p->edges; feature_node; feature_node = feature_node->next) n++;
	for (feature_node = p->faces; feature_node; feature_node = feature_node->next) n++;
	
	return n;
}