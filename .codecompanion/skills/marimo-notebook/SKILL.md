---
name: marimo-notebook
description: Write a marimo notebook in a Python file in the right format.
---

# Notes for marimo Notebooks

## Structure of a python notebook

### Example

marimo uses Python to create notebooks, unlike Jupyter which uses JSON. Here's an example notebook: 

```python
# /// script
# dependencies = [
#     "marimo",
#     "numpy==2.4.3",
# ]
# requires-python = ">=3.14"
# ///

import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


# Define setup cell
with app.setup:
    import numpy as np


@app.cell
def _():
    import marimo as mo

    return mo, np

# Define function cell
@app.function
def calculate_statistics(data):
    """Calculate basic statistics for a dataset"""
    return {
        "mean": np.mean(data),
        "median": np.median(data),
        "std": np.std(data)
    }


@app.cell
def _():
    print("hello world")
    return


@app.cell
def _(np, slider):
    np.array([1,2,3]) + slider.value
    return


@app.cell
def _(mo):
    slider = mo.ui.slider(1, 10, 1, label="number to add")
    slider
    return (slider,)


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
```

Notice how the notebook is structured with functions can represent cell contents. Each cell is defined with the `@app.cell` decorator and the inputs/outputs of the function are the inputs/outputs of the cell.
The decorator `with app.setup` also defines a special kind of cell called setup_cell (which is unique within a notebook) marimo usually takes care of the dependencies between cells automatically. 


### Top level imports

You can import top-level functions and classes defined in a marimo notebook into other Python scripts or notebooks using normal Python syntax, as long as your definitions satisfy the simple criteria described on this page. This makes your notebook code reusable, testable, and easier to edit in text editors of your choice.

For a function or class to be saved at the top level of the notebook file, it must meet the following criteria:

The cell must define just a single function or class.
The defined function or class can only refer to symbols defined in the setup cell, or to other top-level symbols.

In the introductory example, the setup cell is represented as a context manager `app.setup` and the cell that contains `calculate_statistics` is represented as a function decorator `@app.function`. You can now import `calculate_statistics` from other Python scripts or notebooks. There can be no more than one setup cell per notebook.

```python
# In another_script.py
from my_notebook import calculate_statistics

data = [1, 2, 3, 4, 5]
stats = calculate_statistics(data)
print(stats)
```


## State and Reactivity

### Instructions

Variables between cells define the reactivity of the notebook for 99% of the use-cases out there. No special state management needed. Don't mutate objects across cells (e.g., `my_list.append()`); create new objects instead. Avoid `mo.state()` unless you need bidirectional UI sync or accumulated callback state. See [STATE.md](references/STATE.md) for details.

### The DAG

marimo statically analyzes each cell to build a directed acyclic graph:

- **References** = global variables the cell reads (function parameters)
- **Definitions** = global variables the cell creates (return tuple)

When a cell runs, all cells that reference its definitions automatically run. Execution order is determined by the DAG, not cell position on the page.

### Variable Uniqueness

Every global variable must be defined by exactly one cell. This prevents ambiguity in the dependency graph.

If you need the same name in multiple cells, use underscore-prefixed cell-local variables:

```python
@app.cell
def _():
    _temp = expensive_computation()
    result_a = summarize(_temp)
    return (result_a,)

@app.cell
def _():
    _temp = different_computation()  # no conflict, _temp is cell-local
    result_b = summarize(_temp)
    return (result_b,)
```

### Mutations Are Not Tracked

marimo does **not** detect mutations like `.append()`, attribute assignment, or in-place DataFrame operations across cells.

```python
# BAD — mutation in another cell, marimo won't re-run dependents
# Cell 1
items = [1, 2, 3]
# Cell 2
items.append(4)  # invisible to the DAG

# GOOD — create a new variable
# Cell 2
extended = items + [4]
```

Mutations within the same cell that defines the variable are fine:

```python
@app.cell
def _(pd):
    df = pd.DataFrame({"a": [1, 2]})
    df["b"] = [3, 4]  # same cell, fine
    return (df,)
```

### Deleting Cells

Deleting a cell removes its global variables from memory. Cells that referenced those variables become invalidated.


## Key Principle: Keep It Simple

**Show all UI elements always.** Only change the data source in script mode.

- Sliders, buttons, widgets should always be created and displayed
- Don't use try/except for normal control flow

### Good Pattern

```python
# Always show the widget
@app.cell
def _(ScatterWidget, mo):
    scatter_widget = mo.ui.anywidget(ScatterWidget())
    scatter_widget
    return (scatter_widget,)

# Only change data source based on mode
@app.cell
def _(is_script_mode, make_moons, scatter_widget, np, torch):
    if is_script_mode:
        # Use synthetic data for testing
        X, y = make_moons(n_samples=200, noise=0.2)
        X_data = torch.tensor(X, dtype=torch.float32)
        y_data = torch.tensor(y)
        data_error = None
    else:
        # Use widget data in interactive mode
        X, y = scatter_widget.widget.data_as_X_y
        # ... process data ...
    return X_data, y_data, data_error

# Always show sliders - use their .value in both modes
@app.cell
def _(mo):
    lr_slider = mo.ui.slider(start=0.001, stop=0.1, value=0.01)
    lr_slider
    return (lr_slider,)

# Auto-run in script mode, wait for button in interactive
@app.cell
def _(is_script_mode, train_button, lr_slider, run_training, X_data, y_data):
    if is_script_mode:
        # Auto-run with slider defaults
        results = run_training(X_data, y_data, lr=lr_slider.value)
    else:
        # Wait for button click
        if train_button.value:
            results = run_training(X_data, y_data, lr=lr_slider.value)
    return (results,)
```

## Don't Guard Cells with `if` Statements

Marimo's reactivity means cells only run when their dependencies are ready. Don't add unnecessary guards:

```python
# BAD - the if statement prevents the chart from showing
@app.cell
def _(plt, training_results):
    if training_results:  # WRONG - don't do this
        fig, ax = plt.subplots()
        ax.plot(training_results['losses'])
        fig
    return

# GOOD - let marimo handle the dependency
@app.cell
def _(plt, training_results):
    fig, ax = plt.subplots()
    ax.plot(training_results['losses'])
    fig
    return
```

The cell won't run until `training_results` has a value anyway.

## Don't Use try/except for Control Flow

Don't wrap code in try/except blocks unless you're handling a specific, expected exception. Let errors surface naturally.

```python
# BAD - hiding errors behind try/except
@app.cell
def _(scatter_widget, np, torch):
    try:
        X, y = scatter_widget.widget.data_as_X_y
        X = np.array(X, dtype=np.float32)
        # ...
    except Exception as e:
        return None, None, f"Error: {e}"

# GOOD - let it fail if something is wrong
@app.cell
def _(scatter_widget, np, torch):
    X, y = scatter_widget.widget.data_as_X_y
    X = np.array(X, dtype=np.float32)
    # ...
```

Only use try/except when:
- You're handling a specific, known exception type
- The exception is expected in normal operation (e.g., file not found)
- You have a meaningful recovery action

## Cell Output Rendering

Marimo only renders the **final expression** of a cell. Indented or conditional expressions won't render:

```python
# BAD - indented expression won't render
@app.cell
def _(mo, condition):
    if condition:
        mo.md("This won't show!")  # WRONG - indented
    return

# GOOD - final expression renders
@app.cell
def _(mo, condition):
    result = mo.md("Shown!") if condition else mo.md("Also shown!")
    result  # This renders because it's the final expression
    return
```

## Marimo api docs

- For UI elements in marimo [UI.md](references/UI.md) gives you the UI elements available in marimo
- to get more information about a UI element, or 
if the user wants you to use a marimo functionnality, you can use the context7 mcp server to get relevant information from the API documentation

