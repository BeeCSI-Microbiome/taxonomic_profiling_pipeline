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
    # Base string to which filter targets will be appended
    F_STRING="unclassified\|cellular organism"
    if config['filter_targets']:
        F_STRING = '{0}\|{1}'.format(F_STRING, "\|".join(config["filter_targets"]))
    return F_STRING
