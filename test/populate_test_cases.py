import git
import os
import random
import shutil
import string

# Local imports
from filenames import PHASE_SHIFT_INPUT_FILE, PHASE_SHIFT_OUTPUT_FILE

TEST_CASES_BASE_DIR = "./tests/test_cases"


def main():
    current_repo = git.Repo(".")
    original_head = current_repo.active_branch.name

    try:
        current_heads = current_repo.heads.copy()
        for head in current_heads:
            name = head.name
            if name.endswith("-run"):
                case_id = get_case_id_from_name(name)
                case_dir = create_case_dir(case_id)
                inputs_dir = create_inputs_dir(case_dir)
                outputs_dir = create_outputs_dir(case_dir)

                current_repo.heads[name].checkout()
                print(f"Checked out {name}")

                # Create temporary head at penultimate commit
                new_head_name = create_new_head_name()
                assert new_head_name not in [h.name for h in current_repo.heads]
                current_repo.create_head(new_head_name, "HEAD~1")

                # Copy original input from penultimate commit
                current_repo.heads[new_head_name].checkout()
                shutil.copy(PHASE_SHIFT_INPUT_FILE, inputs_dir)

                # Copy expected outputs from latest commit
                current_repo.heads[name].checkout()
                shutil.copy(PHASE_SHIFT_INPUT_FILE, outputs_dir)
                shutil.copy(PHASE_SHIFT_OUTPUT_FILE, outputs_dir)

                # Delete temporary head
                current_repo.delete_head(new_head_name)
    finally:
        current_repo.heads[original_head].checkout()


def get_case_id_from_name(name):
    case_id = name.split("-")[0]
    return case_id


def create_case_dir(case_id):
    assert os.path.isdir(TEST_CASES_BASE_DIR)
    case_dir = "case_" + case_id
    path = create_sub_dir(TEST_CASES_BASE_DIR, case_dir)
    return path


def create_inputs_dir(case_dir):
    path = create_sub_dir(case_dir, "inputs")
    return path


def create_outputs_dir(case_dir):
    path = create_sub_dir(case_dir, "outputs")
    return path


def create_sub_dir(parent, child):
    path = os.path.join(parent, child)
    os.mkdir(path)
    return path


def create_new_head_name():
    alphabet = list(set(string.hexdigits.lower()))
    name = "".join(random.choices(alphabet, k=10))
    return name


if __name__ == "__main__":
    main()
