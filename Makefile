all: build

build:
	cargo build

update:
	cargo update

clean:
	cargo clean

.PHONY: \
	all \
	build \
	clean \
	update
