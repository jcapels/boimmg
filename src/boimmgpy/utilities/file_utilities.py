
def read_conf_file(path):
    """
    reading the configuration file
    It must be assembled in the following manner: key = configuration value
    It will return a dictionary in the following format: {workers : 3, limit : 100}
    """
    res = {}

    with open(path) as file:
        lines = file.readlines()
        for line in lines:
            if '=' in line:
                prop = line.strip()
                l_str = prop.split('=')
                if len(l_str) >2:
                    l_str[1] = "=".join(l_str[1:])

                res[l_str[0].strip()] = l_str[1].strip()

    return res


def change_conf_file(path, dic):
    new_lines = []
    with open(path, "r") as file:
        lines = file.readlines()
        for line in lines:
            if '=' in line:
                prop = line.strip()
                l_str = prop.split('=')
                new_prop = dic[l_str[0].strip()]
                new_line = l_str[0] + "=" + new_prop + "\n"
                new_lines.append(new_line)
            else:
                new_lines.append(line)

    with open(path, "w") as file:
        file.writelines(new_lines)
