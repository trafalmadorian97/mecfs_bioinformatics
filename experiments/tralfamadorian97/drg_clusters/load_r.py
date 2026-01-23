import pyreadr

def go():
    result = pyreadr.read_r("human_meta_final_cluster.Rdata")
    import pdb; pdb.set_trace()
    print(result.keys())

if __name__ == '__main__':
    go()