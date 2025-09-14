# Makefile

# author = "Alejandro Gonzales-Irribarren"
# email = "alejandrxgzi@gmail.com"
# github = "https://github.com/alejandrogzi"
# version: 0.0.1

.PHONY: configure test

configure:
	@bash isopipe/assets/sh/cfg.sh
	@$(MAKE) test

test:
	@echo "INFO: Testing installation..."
	@isopipe/target/release/isopipe --help >/dev/null && echo 'INFO: isopipe installed!' || (echo "ERROR: isopipe run failed!" && exit 1)
	@../isotools/isotools/target/release/isotools --help >/dev/null && echo 'INFO: isopipe installed!' || (echo "ERROR: isotools build failed!" && exit 1)
	@isopipe/assets/rust/orf/target/release/orf --help >/dev/null && echo 'INFO: orf installed!' || (echo "ERROR: orf build failed!" && exit 1)
	@test -f isopipe/assets/rust/orf/tai/.venv/bin/python || (echo "ERROR: Python venv missing!" && exit 1)
	@echo "INFO: All tests passed! You can run the pipeline with: isopipe --help"
