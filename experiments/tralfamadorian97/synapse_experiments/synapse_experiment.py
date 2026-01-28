
import synapseclient
import synapseutils
from synapseclient.models import File, Folder, Project



def go():
    syn = synapseclient.login()
    # project = Project(name="UKB-PPP").get()
    # project.sync_from_synapse(download_file=False)

    # walk returns a generator yielding (dirpath, dirnames, filenames)
    project_or_folder_id = "syn51365303"
    walked_path = synapseutils.walk(syn, project_or_folder_id)

    for dirpath, dirnames, filenames in walked_path:
        for filename in filenames:
            # filename[0] is the name, filename[1] is the Synapse ID
            print(f"File: {filename[0]} (ID: {filename[1]})")

    # dir_mapping = project.map_directory_to_all_contained_files("./")
    # for directory_name, file_entities in dir_mapping.items():
    #     print(f"Directory: {directory_name}")
    #     for file_entity in file_entities:
    #         print(f"\tFile: {file_entity.name}, ID: {file_entity.id}")
    #

def try_download():
    syn = synapseclient.login()
    files = synapseutils.syncFromSynapse(syn, 'syn52363617', path="dummy")
    import pdb; pdb.set_trace()
    synapseutils.syncFromSynapse(syn, 'syn52363617', path="dummy")
    print("yo")

if __name__ == "__main__":
    go()
    # try_download()