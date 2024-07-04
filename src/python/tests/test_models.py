import os
import pytest
import sys
import yt

from numpy.testing import assert_allclose

from pygrackle.utilities.model_tests import model_sets
from pygrackle.utilities.testing import \
    dirname, \
    run_command, \
    temporary_directory

grackle_install_dir = dirname(os.path.abspath(__file__), level=4)
grackle_python_dir = os.path.join(grackle_install_dir, "src", "python")
python_example_dir = os.path.join(grackle_python_dir, "examples")
test_answers_dir = os.path.join(grackle_python_dir, "tests", "example_answers")

generate_results = \
  int(os.environ.get("GENERATE_MODEL_TEST_RESULTS", 0)) == 1

all_sets = []
for model_name, model in model_sets.items():
    for par_index in range(len(model["parameter_sets"])):
        for input_index in range(len(model["input_sets"])):
            all_sets.append((model_name, par_index, input_index))

@pytest.mark.parametrize("model_name, par_index, input_index", all_sets)
def test_model(model_name, par_index, input_index):
    script_path = os.path.join(python_example_dir, f"{model_name}.py")
    command = f"{sys.executable} {script_path} {par_index} {input_index}"
    with temporary_directory():
        rval = run_command(command, timeout=60)
        assert rval

        output_file = f"{model_name}_{par_index}_{input_index}.h5"
        answer_path = os.path.join(test_answers_dir, output_file)

        if generate_results:
            os.rename(output_file, answer_path)
        else:
            assert os.path.exists(answer_path)

            ds1 = yt.load(output_file)
            ds2 = yt.load(answer_path)
            assert ds1.field_list == ds2.field_list

            for field in ds1.field_list:
                err_msg = f"Model mismatch: {model_name}, {par_index}, " + \
                  f"{input_index}: {field}."
                assert_allclose(ds1.data[field], ds2.data[field],
                                atol=0, rtol=1e-8, err_msg=err_msg)
