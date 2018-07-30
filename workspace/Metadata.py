from airtable import Airtable

from . import Workspace


def get_samples(metadata_filter):

    if not is_data_loaded:
        load_data_from_airtable()

    global FASTQ_files
    global samples
    global experiments

    current_experiment = None

    for experiment in experiments:
        if experiment["fields"]["Name"] == Workspace.experiment_name:
            current_experiment = experiment
            break

    if current_experiment is None:
        raise EnvironmentError("Current experiment not set")

    sample_names = []

    for sample in samples:

        if sample["fields"]["Experiment"][0] != current_experiment["id"]:
            continue

        for key, value in metadata_filter.items():

            if key not in sample["fields"]:
                break
            if sample["fields"][key] == value:
                sample_names.append(sample["fields"]["Name"])

    return sample_names


def load_data_from_airtable():

    FASTQ_files_table = Airtable("app4AhX7PjZCeg1ws", "FASTQ Files")
    samples_table = Airtable("app4AhX7PjZCeg1ws", "Samples")
    experiments_table = Airtable("app4AhX7PjZCeg1ws", "Experiments")

    global FASTQ_files
    global samples
    global experiments

    FASTQ_files = FASTQ_files_table.get_all()
    samples = samples_table.get_all()
    experiments = experiments_table.get_all()

    global is_data_loaded

    is_data_loaded = True


is_data_loaded = False
samples = None
experiments = None
FASTQ_files = None