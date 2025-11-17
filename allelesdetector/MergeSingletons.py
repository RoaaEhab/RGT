from typing import List, Tuple, Dict, Any
from rapidfuzz.distance import Levenshtein

def merge_singletons(sorted_geno_list: List[Tuple[str, List[Any]]], top_n: int = 2, max_distance: int =2 , singleton_threshold: int = 1) -> Tuple[List[Tuple[str, List[Any]]], Dict[str, List[Tuple[str, int]]]]:
    
    dominants = sorted_geno_list[:top_n]
    dominant_info = [(name, data, data[5]) for name, data in dominants] 
    mapping: Dict[str, List[Tuple[str, int]]] = {name: [] for name, _ in dominants}
    dominant_abundance_updates: Dict[str, int] = {name: 0 for name, _ in dominants}
    
    for name, data in sorted_geno_list:
        abundance = data[0]
        
        if abundance > singleton_threshold:
            continue 
        raw_dna = data[5]
        best_match_name = None
        min_dist = max_distance + 1
        
        for dom_name, dom_data, dom_raw_dna in dominant_info:
            dist = Levenshtein.distance(raw_dna, dom_raw_dna)
            if dist < min_dist:
                min_dist = dist
                best_match_name = dom_name
        
        if best_match_name and min_dist <= max_distance:
            mapping[best_match_name].append((name, abundance))
            dominant_abundance_updates[best_match_name] += abundance
    
    cleaned: List[Tuple[str, List[Any]]] = []
    for name, data in sorted_geno_list:
        if name in dominant_abundance_updates:
            new_data = data.copy()
            new_data[0] = data[0] + dominant_abundance_updates[name]
            cleaned.append((name, new_data))
        else:
            cleaned.append((name, data.copy()))

    return cleaned, mapping