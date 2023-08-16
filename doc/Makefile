ENV = .venv

all: build

serve:
	${ENV}/bin/mkdocs serve

build: .venv
	${ENV}/bin/mkdocs build

env: .venv

.venv:
	python3 -m venv ${ENV}
	.venv/bin/python3 -m pip install mkdocs
	.venv/bin/python3 -m pip install mkdocs-material
	.venv/bin/python3 -m pip install mkdoxy

clean-env:
	rm -rf ${ENV}

clean-site:
	rm -rf site

clean: clean-env clean-site
