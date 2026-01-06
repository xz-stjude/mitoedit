# MitoEdit Tests

This directory contains pytest tests for the MitoEdit package.

## Running Tests

### Prerequisites
Make sure you're in the MitoEdit root directory:
```bash
cd /path/to/MitoEdit
```

### Install test dependencies
```bash
pip install -e .
```

### Run all tests
```bash
pytest tests/
```

### Run specific test file
```bash
pytest tests/test_findTAL.py -v
```

### Run with coverage
```bash
pytest tests/ --cov=mitoedit --cov-report=html
```

### Run with verbose output
```bash
pytest tests/ -v
```

### Alternative: Run from tests directory
```bash
cd tests/
python -m pytest test_findTAL.py -v
```

## Test Structure

- `conftest.py` - Shared pytest fixtures
- `test_findTAL.py` - Tests for the findTAL module specifically
- `test_talent_tools.py` - Integration tests for the talent_tools module
- `README.md` - This file

## Modern Python Package Test Standards

This test structure follows modern Python packaging standards:

```
MitoEdit/
├── mitoedit/           # Source code package
│   ├── talent_tools/   # Module being tested
│   └── ...
├── tests/              # Test directory (same level as source)
│   ├── __init__.py     # Makes tests a package
│   ├── conftest.py     # Shared fixtures
│   ├── test_*.py       # Test files
│   └── README.md       # Test documentation
├── pyproject.toml      # Package configuration
└── pytest.ini         # Pytest configuration
```

**Key principles:**
- Tests are **outside** the source package
- Tests mirror the source structure
- Each module has corresponding test files
- Shared fixtures in `conftest.py`
- Tests can be run independently of the package

## Test Categories

### Unit Tests
- Test individual functions in isolation
- Use mocks for external dependencies
- Fast execution

### Integration Tests
- Test complete workflows
- Use real files and data
- May be slower

## Writing New Tests

1. Create test files with `test_` prefix
2. Use descriptive test function names starting with `test_`
3. Use fixtures from `conftest.py` for common setup
4. Mock external dependencies when appropriate
5. Add docstrings explaining what each test validates
