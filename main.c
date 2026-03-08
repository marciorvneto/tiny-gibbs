#include <stdio.h>
#include <stdlib.h>

//================================
//
//   Arena
//
//================================

typedef struct {
  char *base;
  size_t offset;
  size_t capacity;
} Arena;

Arena arena_create(size_t capacity) {
  Arena a = {0};
  a.capacity = capacity;
  a.base = malloc(capacity);
  return a;
}

void arena_destroy(Arena *a) { free(a->base); }

void *arena_alloc(Arena *a, size_t amount) {
  size_t aligned = (amount + 7) & ~7;
  if (a->offset + aligned > a->capacity) {
    return NULL;
  }
  void *addr = a->base + a->offset;
  a->offset += aligned;
  return addr;
}

int main() {
  //
  return 0;
}
