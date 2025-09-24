import json
def parse_match_results(file_path):
    """
    Parses the PharmCAT match results from a JSON file.

    Args:
        file_path (str): The path to the JSON file containing PharmCAT match results.

    Returns:
        dict: A dictionary containing the parsed match results.
    """
    with open(file_path, 'r') as file:
        data = json.load(file)
        match_list = filter(lambda x:len(x['diplotypes'])>0,data['results']).tolist()

    return data

def parse_report_results(file_path):
    """
    Parses the PharmCAT report results from a JSON file.

    Args:
        file_path (str): The path to the JSON file containing PharmCAT report results.

    Returns:
        dict: A dictionary containing the parsed report results.
    """
    with open(file_path, 'r') as file:
        data = json.load(file)
        report_list = filter(lambda x:len(x['diplotypes'])>0,data['results']).tolist()

    return data