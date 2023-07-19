all: build

serve:
	mkdocs serve

build:
	mkdocs build

env: .venv

.venv:
	python3 -m venv .venv

clean-env:
	rm -rf .venv

clean-site:
	rm -rf site

clean: clean-env clean-site
