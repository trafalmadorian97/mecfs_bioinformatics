import ieugwaspy
def go():
    # result = ieugwaspy.gwasinfo_files("finn-b-M13_ANKYLOSPON")
    result = ieugwaspy.gwasinfo_files("ukb-b-18194")
    print(result)
    import pdb; pdb.set_trace()
    print("testing")

if __name__ == '__main__':
    go()
