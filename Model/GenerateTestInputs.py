#
# FUNCTIONS FOR GENERATING THE INPUTS FOR TESTING
#

def CreateIDPRanges(file_path):
    id_ranges = [(100, 80), (80, 60), (60, 40), (40, 20), (20, 0)]
    result = {}
    def get_range(pident):
        for high, low in id_ranges:
            if high >= pident > low:
                return f"{high}-{low}"
    with open(file_path, 'r') as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            query_id = parts[0]
            target_id = parts[1]
            print(query_id, target_id)
            if query_id == target_id: continue
            pident = float(parts[2])*100
            print(pident)
            range_key = get_range(pident)
            if query_id not in result:
                result[query_id] = {f"{h}-{l}": [] for h, l in id_ranges}
            result[query_id][range_key].append(target_id)
    return result
