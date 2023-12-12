#include "stddef.h"

typedef struct {
  char **array;
  size_t used;
  size_t size;
} ArrayString;

void initArrayString(ArrayString *a, size_t initialSize) {
  a->array = malloc(initialSize * sizeof(char *));
  for (size_t i = 0; i < initialSize; i++) {
    a->array[i] = NULL;
  }
  a->used = 0;
  a->size = initialSize;
}

void insertArrayString(ArrayString *a, char *element) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = realloc(a->array, a->size * sizeof(char *));
  }
  a->array[a->used++] = strdup(element);
}

void freeArrayString(ArrayString *a) {
  for (size_t i = 0; i < a->used; i++) {
    free(a->array[i]);
  }
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}