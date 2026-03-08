OUT_DIR := ./out
TARGET := $(OUT_DIR)/main
CPPFLAGS := -I./vendor/tinyla
LDLIBS := -g -lm

all: $(TARGET) $(TARGET_BROWSER)

$(TARGET): main.c | $(OUT_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ $(LDLIBS)

$(OUT_DIR):
	@mkdir -p $(OUT_DIR)

clean:
	@rm -rf $(OUT_DIR)
