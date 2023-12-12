#include "stddef.h"

typedef struct {
  float *array;
  size_t used;
  size_t size;
} ArrayFloat;

void initArrayFloat(ArrayFloat *a, size_t initialSize) {
  a->array = malloc(initialSize * sizeof(float));
  a->used = 0;
  a->size = initialSize;
}

void insertArrayFloat(ArrayFloat *a, float element) {
  // a->used is the number of used entries, because a->array[a->used++] updates a->used only *after* the array has been accessed.
  // Therefore a->used can go up to a->size 
  if (a->used == a->size) {
    a->size *= 2;
    a->array = realloc(a->array, a->size * sizeof(float));
  }
  a->array[a->used++] = element;
}

void freeArrayFloat(ArrayFloat *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}