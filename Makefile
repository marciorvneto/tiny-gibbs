OUT_DIR := ./out
OUT_BROWSER_DIR := ./docs
TARGET := $(OUT_DIR)/main
CPPFLAGS := -I./vendor/tinyla
LDLIBS := -g -lm
TARGET_BROWSER := $(OUT_BROWSER_DIR)/main.js

EMCC_FLAGS := -O3 -s ALLOW_MEMORY_GROWTH=1 \
  -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' \
  -s EXPORTED_FUNCTIONS='["_main", "_init","_solve","_malloc","_free"]'

all: $(TARGET) $(TARGET_BROWSER)

$(TARGET): main.c | $(OUT_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ $(LDLIBS)

$(TARGET_BROWSER): main.c | $(OUT_BROWSER_DIR)
	emcc $(EMCC_FLAGS) $(CPPFLAGS) $< -o $@ -lm --embed-file db/nasa9_combustion.dat@db/nasa9_combustion.dat
	@cp ./index.html $(OUT_BROWSER_DIR)/index.html

$(OUT_DIR):
	@mkdir -p $(OUT_DIR)
$(OUT_BROWSER_DIR):
	mkdir -p $(OUT_BROWSER_DIR)

clean:
	@rm -rf $(OUT_DIR)
	@rm -rf $(OUT_BROWSER_DIR)
