import pathlib
import subprocess
import pytest
import miniramses as ram
import numpy as np


@pytest.fixture()
def system_tester():
    class SystemTestRunner:
        def __init__(self):
            # define some useful paths
            self.repo_root = pathlib.Path(__file__).resolve().parents[1]
            self.bin_dir = self.repo_root / "bin"
            self.ranks = 22
            self.system_test_allowed_error = 1.0e-15

        def build_run_and_test(self, build_options, test_name, output_num=2):
            # Build Ramses
            self.build_ramses(build_options)

            # Run Ramses, pass name of namelist file
            self.run_ramses(test_name)

            # Check the results
            self.check_results(test_name, output_num)

        def build_ramses(self, build_options: list[str]):
            # Define make commands and ramses executable name
            make_clean = ["make", "clean"]
            build_command = ["make"] + build_options + ["COMPILER=NVHPC", "MPI=1"]

            dims = next((s for s in build_command if "NDIM" in s))[
                -1
            ]  # find the value of the NDIM option
            ramses_name = f"ramses{dims}d"

            # Create directory for executable and test results. The directory name is the `build_command` command without the leading `make` and with all `=` removed. If the directory already exists then this build has already been done so just return
            self.build_dir = self.bin_dir / ("_".join(build_command[1:])).replace(
                "=", ""
            )
            if self.build_dir.is_dir():
                self.ramses_path = self.build_dir / ramses_name
                return
            self.build_dir.mkdir()

            # make directory
            with open(self.build_dir / "compile.log", "w") as log_file:
                subprocess.run(
                    make_clean,
                    cwd=self.bin_dir,
                    stderr=subprocess.STDOUT,
                    stdout=log_file,
                    check=True,
                )
                subprocess.run(
                    build_command,
                    cwd=self.bin_dir,
                    stderr=subprocess.STDOUT,
                    stdout=log_file,
                    check=True,
                )
                subprocess.run(
                    make_clean,
                    cwd=self.bin_dir,
                    stderr=subprocess.STDOUT,
                    stdout=log_file,
                    check=True,
                )

            # copy ramses exe to build directory
            self.ramses_path = (self.bin_dir / ramses_name).rename(
                self.build_dir / ramses_name
            )

        def run_ramses(self, test_name):  # returns test_dir
            # setup directory for this test
            self.test_dir = self.build_dir / test_name
            self.test_dir.mkdir(exist_ok=True)

            # Run Ramses
            namelist_file = self.repo_root / "tests" / "namelist" / f"{test_name}.nml"
            run_command = [
                "mpirun",
                "-n",
                f"{self.ranks}",
                self.ramses_path.as_posix(),
                namelist_file.as_posix(),
            ]
            with open(self.test_dir / "run.log", "w") as log_file:
                subprocess.run(
                    run_command,
                    cwd=self.test_dir,
                    stderr=subprocess.STDOUT,
                    stdout=log_file,
                    check=True,
                )

        def check_results(self, test_name, output_num=2):
            # Load test & fiducial data
            fiducial_path = (
                self.repo_root / "mini-ramses-test-data" / "system_tests" / test_name
            )
            fiducial_data = ram.rd_cell(output_num, path=fiducial_path.as_posix())
            test_data = ram.rd_cell(output_num, path=self.test_dir.as_posix())

            diff = np.abs(fiducial_data.u - test_data.u)
            l1_error = np.mean(diff, axis=1)
            l2_error = np.linalg.norm(l1_error)

            assert l2_error < self.system_test_allowed_error

    return SystemTestRunner()


def test_system_sedov2d(system_tester):
    build_options = ["NDIM=2", "HYDRO=1"]
    test_name = "sedov2d"

    system_tester.build_run_and_test(build_options, test_name)


def test_system_sedov3d(system_tester):
    build_options = ["NDIM=3", "HYDRO=1"]
    test_name = "sedov3d"

    system_tester.build_run_and_test(build_options, test_name)
