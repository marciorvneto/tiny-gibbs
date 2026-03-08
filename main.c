#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DB_PATH "./db/nasa9_combustion.dat"
#define MAX_LINE_SIZE 1024
#define MAX_ATOMS 128
#define MAX_COMPONENTS 128

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

//================================
//
//   IO
//
//================================

const char *read_db_file(Arena *a, const char *path) {
  FILE *db_file = fopen(path, "r");
  if (db_file == NULL) {
    printf("Could not open the database file\n");
    exit(1);
  }

  if (fseek(db_file, 0, SEEK_END) != 0) {
    fclose(db_file);
    return NULL;
  }

  long size = ftell(db_file);

  rewind(db_file);

  char *buffer = arena_alloc(a, (size_t)size + 1); // Null-terminated
  if (buffer == NULL) {
    printf("Arena ran out of memory when reading database.\n");
    exit(1);
  }

  size_t bytes_read = fread(buffer, 1, size, db_file);
  buffer[bytes_read] = '\0';

  fclose(db_file);

  return buffer;
}

//================================
//
//   DB
//
//================================

typedef struct {
  char symbol[4]; // "C", "H", "O", "N", "Ar"
  int count;
} AtomEntry;

typedef struct {
  char name[16];
  int n_elements;
  double mol_weight;
  double hf298; // J/mol
  double T_low, T_mid, T_high;
  AtomEntry atoms[5]; // max 5 elements (C,H,O,N,Ar)
  double low[9];      // NASA-9 coeffs, low-T range
  double high[9];     // NASA-9 coeffs, high-T range
} DBComponent;

typedef struct {
  DBComponent components[MAX_COMPONENTS];
  size_t num_components;
} DBData;

//================================
//
//   DB Parser
//
//================================

typedef struct {
  const char *string;
  size_t pointer;
  size_t line;
} Parser;

typedef struct {
  const char *start;
  size_t length;
  size_t line;
} LineView;

Parser create_parser(const char *string) {
  Parser p = {0};
  p.line = 1;
  p.string = string;
  return p;
}

// Returns the number of bytes read
int read_line(Parser *p, LineView *line) {
  if (p->string[p->pointer] == '\0')
    return 0;

  size_t start = p->pointer;
  size_t end = start;

  while (p->string[end] != '\n' && p->string[end] != '\0') {
    end++;
  }

  line->start = p->string + start;
  line->length = end - start;
  line->line = p->line;

  p->pointer = end; // Default
  if (p->string[end] == '\n') {
    p->line++;
    p->pointer = end + 1;
  }
  if (p->string[end] == '\0') {
    return 0;
  }

  return end - start;
}

void parse_component_atoms(const char *line_str, DBComponent *c) {

  const char *cursor = line_str;
  int n;
  for (size_t i = 0; i < c->n_elements; i++) {
    sscanf(cursor, "%3s %d%n", c->atoms[i].symbol, &c->atoms[i].count, &n);
    cursor += n;
  }
}

void parse_db(Parser *p, DBData *db) {
  LineView line = {0};
  char line_str[MAX_LINE_SIZE];
  int num_components = -1;

  // Skip comments and get to the number of components
  while (p->string[p->pointer] != '\0') {
    read_line(p, &line);
    if (line.start[0] == '#')
      continue;
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    num_components = atoi(line_str);
    break;
  }

  for (int i = 0; i < num_components; i++) {
    DBComponent *c = &db->components[db->num_components];

    // Line 1: name  n_elements  mol_weight  hf298  T_low  T_mid  T_high
    read_line(p, &line);
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    sscanf(line_str, "%15s %d %lf %lf %lf %lf %lf", c->name, &c->n_elements,
           &c->mol_weight, &c->hf298, &c->T_low, &c->T_mid, &c->T_high);

    // Line 2: element pairs
    read_line(p, &line);
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    parse_component_atoms(line_str, c);

    // Line 3: low-T coefficients
    read_line(p, &line);
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    sscanf(line_str, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &c->low[0],
           &c->low[1], &c->low[2], &c->low[3], &c->low[4], &c->low[5],
           &c->low[6], &c->low[7], &c->low[8]);

    // Line 4: high-T coefficients
    read_line(p, &line);
    memcpy(line_str, line.start, line.length);
    line_str[line.length] = '\0';
    sscanf(line_str, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &c->high[0],
           &c->high[1], &c->high[2], &c->high[3], &c->high[4], &c->high[5],
           &c->high[6], &c->high[7], &c->high[8]);

    db->num_components++;
  }
}

//================================
//
//   Thermodynamics
//
//================================

int main() {
  Arena a = arena_create(1024 * 1024); // 1MB

  const char *contents = read_db_file(&a, DB_PATH);
  Parser p = create_parser(contents);

  DBData db = {0};
  parse_db(&p, &db);

  arena_destroy(&a);
  return 0;
}
