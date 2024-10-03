import warnings
import functools # necessary to Sphinx correctly autodoc

def get_runtime_warning_act_decorator(runtime_warning_action="ignore"):
	"""
	Returns the decorator wraping a function with given reaction to RuntimeWarning

	Parameters:
	-----------
	runtime_warning_action : str
		What should be done with RuntimeWarning, often throwing by optimizers.
		Possible values:
		default         # Show all warnings (even those ignored by default)
		ignore          # Ignore all warnings
		error           # Convert all warnings to errors

	Returns:
	--------
	runtime_warning_act : function
		A decorator function that applies the desired warning handling to a wrapped function.
	"""
	def runtime_warning_act(f):
		@functools.wraps(f) # necessary to Sphinx correctly autodoc
		def wrapper(*args, **kwargs):
			with warnings.catch_warnings():
				warnings.filterwarnings(runtime_warning_action, category=RuntimeWarning)
				return f(*args, **kwargs)
		return wrapper
	return runtime_warning_act