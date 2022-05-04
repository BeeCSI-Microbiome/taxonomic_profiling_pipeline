# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#             Functions             #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

def retain(flag, path):
    """Returns given path if flag is true, else returns temp(path)
       which causes the output to be deleted after it is no longer needed"""
    if flag: return path
    else: return temp(path)

def get_filter_string():
    """Creates the filter string based on targets specified in config file"""
    # Get filter target taxa
    FILTER_TARGETS = config['filter_targets']
    BASE_F_STRING="unclassified\|cellular organism"
    if FILTER_TARGETS:
        print("Found filter targets in the config file")
        F_STRING = BASE_F_STRING + "\|" + "\|".join(FILTER_TARGETS)
    print("Filtering with the following filter string:\n\t'{}'".format(F_STRING))
    return F_STRING