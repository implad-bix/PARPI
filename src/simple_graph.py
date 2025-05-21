import os,sys
from collections import defaultdict


def parse_gfa(file_path):
    sequences = {}
    edges = []
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if fields[0] == 'S':
                seq_name = fields[1]
                sequence = fields[2]
                other_attrs = []
                depth = None
                for i in fields[3:]:
                    if i.startswith('DP:f:') or i.startswith('dp:f:'):
                        depth = round(float(i[5:]), 1)
                    else:
                        other_attrs.append(i)
                sequences[seq_name] = {
                    'sequence': sequence,
                    'depth': depth,
                    'other_attributes': other_attrs
                }
            elif fields[0] == 'L':
                from_node = fields[1]
                from_node_orientation = fields[2]
                to_node = fields[3]
                to_node_orientation = fields[4]
                match_info = fields[5]
                edges.append({
                    'from_node': from_node,
                    'from_node_orientation': from_node_orientation,
                    'to_node': to_node,
                    'to_node_orientation': to_node_orientation,
                    'match_info': match_info
                })
    return sequences, edges


def remove_redundant_nodes(edges):
    edges = [edge for edge in edges if not (edge['from_node'] == edge['to_node'] and
                                            edge['from_node_orientation'] != edge['to_node_orientation'])]
    return edges


def find_self_connected_nodes(edges):
    self_connected_nodes = set()
    for edge in edges:
        if edge['from_node'] == edge['to_node'] and edge['from_node_orientation'] == edge['to_node_orientation']:
            self_connected_nodes.add(edge['from_node'])
    return self_connected_nodes


def find_multijointed_nodes(sequences, edges):
    # Use dictionaries to store the left and right connections for each node
    left_connections = {}
    right_connections = {}

    # Iterate through all the edges to count left and right connections for each node
    for edge in edges:
        from_node = edge['from_node']
        from_node_orientation = edge['from_node_orientation']
        to_node = edge['to_node']
        to_node_orientation = edge['to_node_orientation']
        # If a node is in from_node, its right connection is from_node_orientation = +; left connection is
        # from_node_orientation = -
        if from_node_orientation == '-':
            if from_node in list(left_connections.keys()) :
                left_connections[from_node].append(to_node)
            else:
                left_connections[from_node] = []
                left_connections[from_node].append(to_node)
        elif from_node_orientation == '+':
            if from_node in list(right_connections.keys()) :
                right_connections[from_node].append(to_node)
            else:
                right_connections[from_node] = []
                right_connections[from_node].append(to_node)
        if to_node_orientation == '-':
            if to_node in list(right_connections.keys()):
                right_connections[to_node].append(from_node)
            else:
                right_connections[to_node] = []
                right_connections[to_node].append(from_node)
        elif to_node_orientation == '+':
            if to_node in list(left_connections.keys()):
                left_connections[to_node].append(from_node)
            else:
                left_connections[to_node] = []
                left_connections[to_node].append(from_node)
    # Find multijointed nodes and double-loop nodes
    all_multijointed_nodes = set()
    two_nodes_loop = set()
    self_loop_nodes = set()
    hover_nodes = set()
    two_nodes_loop_mapping = defaultdict(set)  # Use a defaultdict to store sets of common nodes
    node_to_delete = []

    for node in list(sequences.keys()):
        # Count the number of left and right connections
        left_count = len(left_connections.get(node, []))
        right_count = len(right_connections.get(node, []))

        # Check if it's a multijointed node (both left and right have 2 or more connections)
        if left_count >= 2 and right_count >= 2:
            # Debugging: Output the connection counts for each node
            all_multijointed_nodes.add(node)

            # Check if it's a double-loop node: any common connections between left and right
            left_neighbors = set(left_connections.get(node, []))
            right_neighbors = set(right_connections.get(node, []))
            common_neighbors = left_neighbors & right_neighbors  # Find the intersection of left and right
            if len(left_neighbors) == 2 or len(right_neighbors) == 2:
                if list(left_neighbors)[0] == list(left_neighbors)[1] or list(right_neighbors)[0] == list(
                        right_neighbors)[1]:
                    hover_nodes.add(node)
            if common_neighbors and len(left_neighbors) == 2 and len(right_neighbors) == 2:  # If there are common
                # neighbors, it's a double-loop node
                if len(common_neighbors) == 1:
                    # If only one common neighbor and it's the node itself, add it to self_connected_nodes
                    if list(common_neighbors)[0] == node:
                        self_loop_nodes.add(node)
                    else:
                        # If only one common neighbor and it's not the node itself, it's a double-loop node
                        two_nodes_loop.add(node)
                        two_nodes_loop_mapping[node].update(common_neighbors)  # Add common neighbor to the set
                '''elif len(common_neighbors) == 2:
                    # If multiple common neighbors exist
                    for neighbor in list(common_neighbors):
                        if neighbor == node:
                            node_to_delete.append(node)
                            #edge_to_delete.append(node)
                        else:
                            # Otherwise, add the neighbor to double_loop_mapping
                            two_nodes_loop_mapping[node].add(neighbor)
                    # Always add node to double_loop_nodes for multiple common neighbors case
                    two_nodes_loop.add(node)'''

    # Ordinary multijointed nodes are all multijointed nodes minus the double-loop nodes
    multijointed_nodes = all_multijointed_nodes - two_nodes_loop - self_loop_nodes - hover_nodes
    edge_to_delete = []
    if node_to_delete:
        for node in node_to_delete:
            for edge in edges:
                if edge['from_node'] == edge['to_node'] == node and edge['from_node_orientation'] == edge['to_node_orientation']:
                    edge_to_delete.append(edge)
    if edge_to_delete:
        for edge in reversed(edge_to_delete):
            if edge in edges:
                edges.remove(edge)
    return multijointed_nodes, two_nodes_loop, two_nodes_loop_mapping, left_connections, right_connections, edges, self_loop_nodes


def solve_self_loop_nodes(self_loop_nodes, sequences, edges, left_connections, right_connections):
    for node in self_loop_nodes:
        result = None
        node_depth = sequences[node]['depth']
        left_node = next((i for i in left_connections[node] if i != node), None)
        right_node = next((i for i in right_connections[node] if i != node), None)
        if left_node is not None and right_node is not None:
            left_node_depth = sequences[left_node]['depth']
            right_node_depth = sequences[right_node]['depth']
            avg_depth = (left_node_depth + right_node_depth) / 2
            result = round(node_depth / avg_depth)
        if result and result > 1:
            original_node = sequences[node]
            original_sequence = original_node['sequence']
            original_other_attrs = original_node['other_attributes']
            original_depth = original_node['depth']
            avg_depth = round(original_depth / result, 1)
            del sequences[node]
            for i in range(1, result + 1):
                copy_name = f"{node}_copy{i}"
                sequences[copy_name] = {
                    'sequence': original_sequence,
                    'depth': avg_depth,
                    'other_attributes': original_other_attrs
                }
            # Assuming edges is the original list of edges, and result is the number of copies
            new_edges = []
            # Traverse the original edges
            for edge in edges:
                from_node = edge['from_node']
                to_node = edge['to_node']
                # Handle the first type of edge: from_node = left_node, to_node = node, or from_node = node, to_node = left_node
                if from_node == left_node and to_node == node:
                    # Replace node with node_copy1
                    new_edges.append({
                        'from_node': from_node,
                        'from_node_orientation': edge['from_node_orientation'],
                        'to_node': f"{node}_copy1",
                        'to_node_orientation': edge['to_node_orientation'],
                        'match_info': edge['match_info']
                    })
                elif from_node == node and to_node == left_node:
                    # Replace node with node_copy1
                    new_edges.append({
                        'from_node': f"{node}_copy1",
                        'from_node_orientation': edge['from_node_orientation'],
                        'to_node': left_node,
                        'to_node_orientation': edge['to_node_orientation'],
                        'match_info': edge['match_info']
                    })
                # Handle the third type of edge: from_node = right_node, to_node = node, or from_node = node, to_node = right_node
                elif from_node == right_node and to_node == node:
                    # Replace node with node_copy{result}
                    new_edges.append({
                        'from_node': right_node,
                        'from_node_orientation': edge['from_node_orientation'],
                        'to_node': f"{node}_copy{result}",
                        'to_node_orientation': edge['to_node_orientation'],
                        'match_info': edge['match_info']
                    })
                elif from_node == node and to_node == right_node:
                    # Replace node with node_copy{result}
                    new_edges.append({
                        'from_node': f"{node}_copy{result}",
                        'from_node_orientation': edge['from_node_orientation'],
                        'to_node': right_node,
                        'to_node_orientation': edge['to_node_orientation'],
                        'match_info': edge['match_info']
                    })
                # Handle the second type of edge: from_node = node, to_node = node (self-loop)
                elif from_node == node and to_node == node:
                    # Create the necessary copy edges based on the result count
                    for i in range(1, result):
                        new_edges.append({
                            'from_node': f"{node}_copy{i}",
                            'from_node_orientation': edge['from_node_orientation'],
                            'to_node': f"{node}_copy{i + 1}",
                            'to_node_orientation': edge['to_node_orientation'],
                            'match_info': edge['match_info']
                        })
                # For all other cases, retain the original edge
                else:
                    pass # new_edges.append(edge)
            # Remove the original 3 edges that need to be replaced
            edges = [edge for edge in edges if
                     edge['from_node'] not in [left_node, node, right_node] or edge['to_node'] not in [left_node, node, right_node] ]
            # Add all the newly generated edges back to the original edges list
            edges.extend(new_edges)
        elif result and result <= 1:
            # Remove edges where from_node == node or to_node == node
            edges = [edge for edge in edges if edge['from_node'] != node and edge['to_node'] != node]
    return sequences, edges


def solve_two_nodes_loop(two_nodes_loop, two_nodes_loop_mapping, sequences, edges, left_connections, right_connections):
    for node in two_nodes_loop:
        #print(node)
        if len(list(two_nodes_loop_mapping[node])) == 1:
            loop_node = list(two_nodes_loop_mapping[node])[0]
            node_depth = sequences[node]['depth']
            #print(node_depth)
            left_node = next((i for i in left_connections[node] if i != loop_node), None)
            right_node = next((i for i in right_connections[node] if i != loop_node), None)
            if left_node is not None and right_node is not None:
                left_node_depth = sequences[left_node]['depth']
                right_node_depth = sequences[right_node]['depth']
                avg_depth = (left_node_depth + right_node_depth) / 2
                #print(left_node_depth)
                #print(right_node_depth)
                result = max(2, round(node_depth / avg_depth))
                original_node = sequences[node]
                original_sequence = original_node['sequence']
                original_other_attrs = original_node['other_attributes']
                original_depth = original_node['depth']
                original_loop_node = sequences[loop_node]
                original_loop_node_sequence = original_loop_node['sequence']
                original_loop_node_other_attrs = original_loop_node['other_attributes']
                original_loop_node_depth = original_loop_node['depth']
                new_node_avg_depth = round(original_depth / result, 1)
                new_loop_node_avg_depth = round(original_loop_node_depth / (result - 1), 1)
                del sequences[node]
                del sequences[loop_node]
                for i in range(1, result + 1):
                    copy_name = f"{node}_copy{i}"
                    sequences[copy_name] = {
                        'sequence': original_sequence,
                        'depth': new_node_avg_depth,
                        'other_attributes': original_other_attrs
                    }
                for i in range(1, result):
                    copy_name = f"{loop_node}_copy{i}"
                    sequences[copy_name] = {
                        'sequence': original_loop_node_sequence,
                        'depth': new_loop_node_avg_depth,
                        'other_attributes': original_loop_node_other_attrs
                    }
                # Assuming edges is the original list of edges, and result is the number of copies
                new_edges = []
                # Traverse the original edges
                past_node_orien = ''
                for edge in edges:
                    from_node = edge['from_node']
                    from_node_orientation = edge['from_node_orientation']
                    to_node = edge['to_node']
                    to_node_orientation = edge['to_node_orientation']
                    match_info = edge['match_info']
                    if from_node == left_node and to_node == node:
                        new_edges.append({
                            'from_node': from_node,
                            'from_node_orientation': from_node_orientation,
                            'to_node': f"{node}_copy1",
                            'to_node_orientation': to_node_orientation,
                            'match_info': match_info
                        })
                    elif from_node == node and to_node == left_node:
                        if to_node_orientation == '+':
                            this_from_orien = '-'
                        else:
                            this_from_orien = '+'
                        if from_node_orientation == '+':
                            this_to_orien = '-'
                        else:
                            this_to_orien = '+'
                        new_edges.append({
                            'from_node': to_node,
                            'from_node_orientation': this_from_orien,
                            'to_node': f"{node}_copy1",
                            'to_node_orientation': this_to_orien,
                            'match_info': match_info
                        })
                    past_node_orien = from_node_orientation
                for edge2 in edges:
                    from_node = edge2['from_node']
                    from_node_orientation = edge2['from_node_orientation']
                    to_node = edge2['to_node']
                    to_node_orientation = edge2['to_node_orientation']
                    match_info = edge2['match_info']
                    if from_node == node and to_node == loop_node and from_node_orientation == past_node_orien:
                        for i in range(1, result):
                            new_edges.append({
                                'from_node': f"{node}_copy{i}",
                                'from_node_orientation': from_node_orientation,
                                'to_node': f"{loop_node}_copy{i}",
                                'to_node_orientation': to_node_orientation,
                                'match_info': match_info
                            })
                            new_edges.append({
                                'from_node': f"{loop_node}_copy{i}",
                                'from_node_orientation': to_node_orientation,
                                'to_node': f"{node}_copy{i + 1}",
                                'to_node_orientation': from_node_orientation,
                                'match_info': match_info
                            })
                    elif from_node == loop_node and to_node == node and to_node_orientation != past_node_orien:
                        if to_node_orientation == '+':
                            second_from_orien = '-'
                        else:
                            second_from_orien = '+'
                        if from_node_orientation == '+':
                            second_to_orien = '-'
                        else:
                            second_to_orien = '+'
                        for i in range(1, result):
                            new_edges.append({
                                'from_node': f"{node}_copy{i}",
                                'from_node_orientation': second_from_orien,
                                'to_node': f"{loop_node}_copy{i}",
                                'to_node_orientation': second_to_orien,
                                'match_info': match_info
                            })
                            new_edges.append({
                                'from_node': f"{loop_node}_copy{i}",
                                'from_node_orientation': second_to_orien,
                                'to_node': f"{node}_copy{i + 1}",
                                'to_node_orientation': second_from_orien,
                                'match_info': match_info
                            })
                for edge3 in edges:
                    from_node = edge3['from_node']
                    from_node_orientation = edge3['from_node_orientation']
                    to_node = edge3['to_node']
                    to_node_orientation = edge3['to_node_orientation']
                    match_info = edge3['match_info']
                    if from_node == node and to_node == right_node and from_node_orientation == past_node_orien:
                        new_edges.append({
                            'from_node': f"{node}_copy{result}",
                            'from_node_orientation': from_node_orientation,
                            'to_node': right_node,
                            'to_node_orientation': to_node_orientation,
                            'match_info': match_info
                        })
                    elif from_node == right_node and to_node == node and to_node_orientation != past_node_orien:
                        if to_node_orientation == '+':
                            third_from_orien = '-'
                        else:
                            third_from_orien = '+'
                        if from_node_orientation == '+':
                            third_to_orien = '-'
                        else:
                            third_to_orien = '+'
                        new_edges.append({
                            'from_node': f"{node}_copy{result}",
                            'from_node_orientation': third_from_orien,
                            'to_node': right_node,
                            'to_node_orientation': third_to_orien,
                            'match_info': match_info
                        })
                # Remove the original 3 edges that need to be replaced
                edges = [edge for edge in edges if edge['from_node'] not in [left_node, node, loop_node, right_node] or
                         edge['to_node'] not in [left_node, node, loop_node, right_node]]
                # Add all the newly generated edges back to the original edges list
                edges.extend(new_edges)
    return sequences, edges


def BLAST_running(makeblastdb_path, BLASTn_path, ref_fasta, out_dir, gfa_file, threads):
    gfa2fasta = os.path.join(out_dir, 'tmp', 'for_blast.fasta')
    os.makedirs(os.path.dirname(gfa2fasta), exist_ok=True)
    with open(gfa_file, 'r') as gfa:
        with open(gfa2fasta, 'a+') as fasta:
            for line in gfa:
                # Skip header lines or empty lines
                if line.startswith('H') or not line.strip():
                    continue
                # Process the 'S' lines that contain sequences
                if line.startswith('S'):
                    # Split the line by tabs
                    parts = line.strip().split('\t')
                    # The first part is the node name (ID)
                    node_id = parts[1]
                    # The second part is the sequence
                    sequence = parts[2]
                    length = len(sequence)
                    # Write the sequence to the FASTA file
                    fasta.write(f">{node_id}#{length}\n{sequence}\n")
    blast_results = os.path.join(out_dir, 'tmp', 'blast_out.txt')
    os.system(f'{makeblastdb_path} -in {ref_fasta} -out {ref_fasta} -dbtype nucl')
    os.system(f'{BLASTn_path} -query {gfa2fasta} -db {ref_fasta} -evalue 1e-10 -out {blast_results} -outfmt 6 '
              f'-num_threads {threads} -max_target_seqs 1000000000 -qcov_hsp_perc 0.8')
    return blast_results


def process_blast_results(blast_file):
    blast_dict = {}
    # Open and read the file
    with open(blast_file, 'r') as file:
        for line in file:
            # Split the line by tab or space
            columns = line.strip().split()

            if len(columns) < 11:
                continue  # Skip lines that don't have enough columns
            # Extract query ID (first column) and blast length (4th column)
            query_id = columns[0].split('#')[0]  # Get the part after the '#'
            blast_length = int(columns[3])  # Length of the BLAST hit
            query_length = int(columns[0].split('#')[1])  # Length of the query

            # Only keep hits where BLAST length is >= 80% of the query length
            if blast_length >= 0.8 * query_length:
                contig_id = columns[1]  # Contig ID (second column)
                start_location = int(columns[8])  # Start position (9th column)
                end_location = int(columns[9])  # End position (10th column)

                # Make sure start_location < end_location
                if start_location > end_location:
                    start_location, end_location = end_location, start_location

                # Update the dictionary with the extracted information
                if query_id in list(blast_dict.keys()):
                    blast_dict[query_id].append((contig_id, start_location, end_location))
                else:
                    blast_dict[query_id] = []
                    blast_dict[query_id].append((contig_id, start_location, end_location))
    return blast_dict


def solve_multijointed_nodes(multijointed_nodes, left_connections, right_connections, sequences, edges,
                                          makeblastdb_path, BLASTn_path, ref_fasta, out_dir, gfa_file, threads):
    blast_results = BLAST_running(makeblastdb_path, BLASTn_path, ref_fasta, out_dir, gfa_file, threads)
    blast_dict = process_blast_results(blast_results)
    for node in multijointed_nodes:
        print('this node is ' + node)
        print('left_connections:')
        print(left_connections[node])
        left_node = left_connections[node]
        right_node = right_connections[node]
        print('right_connections:')
        print(right_connections[node])
        if len(left_node) >= 2 and len(right_node) >= 2:
            copy = min(len(left_node), len(right_node))
            print(copy)
            num = 1
            while copy > 1:
                node_depth = sequences[node]['depth']
                copy_name = f'{node}_copy{num}'
                new_depth = (node_depth) / copy
                new_node = sequences[node]
                new_sequence = new_node['sequence']
                new_other_attrs = new_node['other_attributes']
                sequences[copy_name] = {
                    'sequence': new_sequence,
                    'depth': new_depth,
                    'other_attributes': new_other_attrs
                }
                print(copy_name)
                for item in list(blast_dict[node]):
                    print(item)
                    min_distance_left = float('inf')  # Initialize to a very large number
                    closest_left = None
                    min_distance_right = float('inf')  # Initialize to a very large number
                    closest_right = None
                    this_node_anchor = item[0]
                    this_node_start = item[1]
                    this_node_end = item[2]
                    for left in left_node:
                        print(left)
                        if left in list(blast_dict.keys()):
                            print(blast_dict[left])
                            for i in blast_dict[left]:
                                this_left_anchor = i[0]
                                this_left_start = i[1]
                                this_left_end = i[2]
                                if this_node_anchor == this_left_anchor:
                                    min_distance_candidate = min(abs(this_left_start - this_node_end),
                                                                 abs(this_node_start - this_left_end))
                                    print(min_distance_candidate)
                                    if min_distance_candidate < min_distance_left and min_distance_candidate < 1000:
                                        min_distance_left = min_distance_candidate
                                        closest_left = left
                    print('closest_left:')
                    print(closest_left)
                    for right in right_node:
                        if right in list(blast_dict.keys()):
                            for i in blast_dict[right]:
                                this_right_anchor = i[0]
                                this_right_start = i[1]
                                this_right_end = i[2]
                                if this_node_anchor == this_right_anchor:
                                    min_distance_candidate = min(abs(this_right_start - this_node_end),
                                                                 abs(this_node_start - this_right_end))
                                    if min_distance_candidate < min_distance_right and min_distance_candidate < 1000:
                                        min_distance_right = min_distance_candidate
                                        closest_right = right
                    print('closest_right:')
                    print(closest_right)
                    new_edges = []
                    # modify edges
                    for edge in edges:
                        from_node = edge['from_node']
                        from_node_orientation = edge['from_node_orientation']
                        to_node = edge['to_node']
                        to_node_orientation = edge['to_node_orientation']
                        match_info = edge['match_info']
                        if (from_node == closest_left or from_node == closest_right) and to_node == node:
                            new_edges.append({
                                'from_node': from_node,
                                'from_node_orientation': from_node_orientation,
                                'to_node': f"{node}_copy{num}",
                                'to_node_orientation': to_node_orientation,
                                'match_info': match_info
                            })
                        elif from_node == node and (to_node == closest_left or to_node == closest_right):
                            new_edges.append({
                                'from_node': f"{node}_copy{num}",
                                'from_node_orientation': from_node_orientation,
                                'to_node': to_node,
                                'to_node_orientation': to_node_orientation,
                                'match_info': match_info
                            })
                        # Remove the original 3 edges that need to be replaced
                    edges = [edge for edge in edges if
                                edge['from_node'] not in [closest_left, node, closest_right] or
                                edge['to_node'] not in [closest_left, node, closest_right]]
                    # Add all the newly generated edges back to the original edges list
                    #print(edges)
                    edges.extend(new_edges)
                    #print(edges)
                    if left_node is not None and closest_left in left_node:
                        left_node = [i for i in left_node if i != closest_left]
                    if right_node is not None and closest_right in right_node:
                        right_node = [i for i in right_node if i != closest_right]
                    blast_dict[node] = [i for i in blast_dict[node] if i != item]
                    break
                copy -= 1
                num += 1
                print(copy)
                break
    return sequences, edges


def write_to_gfa(sequences, edges, gfa_filename):
    with open(gfa_filename, 'w') as f:
        # Write the header line (optional but common)
        f.write('H\tVN:1.0\n')

        # Write the sequence lines
        for seq_name, seq_data in sequences.items():
            sequence = seq_data['sequence']
            depth = seq_data['depth']
            other_attrs = seq_data['other_attributes']
            # Adding depth and other attributes as extra information
            f.write(f"S\t{seq_name}\t{sequence}\tDP:f:{depth}")
            for i in other_attrs:
                f.write(f"\t{i}")
            f.write(f"\n")

        # Write the edge lines
        for edge in edges:
            from_node = edge['from_node']
            to_node = edge['to_node']
            from_node_orientation = edge['from_node_orientation']
            to_node_orientation = edge['to_node_orientation']
            match_info = edge['match_info']
            f.write(f"L\t{from_node}\t{from_node_orientation}\t{to_node}\t{to_node_orientation}\t{match_info}\n")


file_path = sys.argv[1]
out_gfa_path = sys.argv[2]
ref_fasta = sys.argv[3]
out_dir = '/home/user001/test/'
threads = 10
makeblastdb_path = '/home/user001/mambaforge/envs/jll/bin//makeblastdb'
BLASTn_path = '/home/user001/mambaforge/envs/jll/bin//blastn'
sequences, edges = parse_gfa(file_path)
edges = remove_redundant_nodes(edges)
multijointed_nodes, two_nodes_loop, two_nodes_loop_mapping, left_connections, right_connections, edges, self_loop_nodes = \
    find_multijointed_nodes(sequences, edges)
# solve_self_loop_nodes
sequences, edges = solve_self_loop_nodes(self_loop_nodes, sequences, edges, left_connections, right_connections)
# solve_two_nodes_loop
sequences, edges = solve_two_nodes_loop(two_nodes_loop, two_nodes_loop_mapping, sequences, edges, left_connections,
                                        right_connections)
# solve_multijointed_nodes
sequences, edges = solve_multijointed_nodes(multijointed_nodes, left_connections, right_connections, sequences, edges,
                                          makeblastdb_path, BLASTn_path, ref_fasta, out_dir, file_path, threads)


write_to_gfa(sequences, edges, out_gfa_path)
'''print("Sequences:")
for seq_name, seq_info in sequences.items():
    print(f"Name: {seq_name}, Sequence: {seq_info['sequence'][:100]}, depth: {seq_info['depth']}, "
          f"Attributes: {seq_info['other_attributes']}")

print("\nEdges:")
for edge in edges:
    print(f"from_node: {edge['from_node']}, from_node_orientation: {edge['from_node_orientation']}, "
          f"to_node: {edge['to_node']}, to_node_orientation: {edge['to_node_orientation']}, match_info: {edge['match_info']}")

# Print the results
print("self_loop_nodes:", self_loop_nodes)
print("multijointed_nodes:", multijointed_nodes)
print("two_nodes_loop:", two_nodes_loop)
print("two_nodes_loop_mapping:", two_nodes_loop_mapping)'''
