#include "paths.h"

// Return the relative path corresponding to the given values:
//     data/(system_name)/(file_type)/(file_name)
// The returned value must be freed using bcstrfree().
char *cwannier_data_path(char *system_name, char *file_type, char *file_prefix, char *file_name) {
    bstring str_path = bfromcstr("data/");
    bcatcstr(str_path, system_name);
    bcatcstr(str_path, "/");
    bcatcstr(str_path, file_type);
    bcatcstr(str_path, "/");
    bcatcstr(str_path, file_prefix);
    bcatcstr(str_path, file_name);

    char *cstr_path = bstr2cstr(str_path, 'N');

    bdestroy(str_path);
    return cstr_path;
}
