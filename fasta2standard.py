#!/usr/bin/env python3
"""
FASTA to HIVDiversity Group Standard Compliance Converter

This script reads a JSON transformation file exported from fastaCompliance.html
and applies those transformations to a FASTA file to convert it to standard compliance.

Usage (single file):
    python fasta2standard.py <json_file> <input_fasta> [output_fasta] [-quiet]

Usage (directory):
    python fasta2standard.py <json_file> <input_directory> [-quiet]

Examples:
    python fasta2standard.py V806_transformations.json input.fasta output_compliant.fasta
    python fasta2standard.py V806_transformations.json ./fasta_files/
    python fasta2standard.py V806_transformations.json ./fasta_files/ -quiet
"""

import json
import re
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython", file=sys.stderr)
    sys.exit(1)

# Import validation functions from fasta_naming
try:
    from py.fasta_naming import parse_fasta_filename, parse_sequence_name, FastaFileName, SequenceName
except ImportError:
    print("Error: Could not import fasta_naming module. Ensure py/fasta_naming.py exists.", file=sys.stderr)
    sys.exit(1)


def parse_fields(text: str, field_delim: str, subfield_delim: str) -> List[Dict[str, any]]:
    """
    Parse text into fields and subfields based on delimiters.
    
    Args:
        text: Text to parse
        field_delim: Field delimiter (e.g., '_')
        subfield_delim: Subfield delimiter (e.g., '-')
        
    Returns:
        List of parts, each with 'text', 'type' ('field', 'subfield', 'delimiter')
    """
    # Strip file extensions
    text = re.sub(r'\.(fasta|fa|fas)$', '', text, flags=re.IGNORECASE)
    
    parts = []
    if not text:
        return parts
    
    # Split by field delimiter
    fields = text.split(field_delim)
    
    for i, field in enumerate(fields):
        if i > 0:
            parts.append({'text': field_delim, 'type': 'delimiter'})
        
        if not field:
            continue
        
        # Split by subfield delimiter
        subfields = field.split(subfield_delim)
        
        for j, subfield in enumerate(subfields):
            if j > 0:
                parts.append({'text': subfield_delim, 'type': 'subdelimiter'})
            
            if subfield:
                parts.append({'text': subfield, 'type': 'subfield' if j > 0 else 'field'})
    
    return parts


def convert_replacement_string(replacement: str) -> str:
    """
    Convert JavaScript-style replacement string to Python format.
    
    JavaScript uses $1, $2, etc. for backreferences, while Python uses \1, \2, etc.
    Also handles literal $ signs (converts $$ to $ for Python, then escapes remaining $).
    
    Args:
        replacement: JavaScript-style replacement string (e.g., "$1" or "prefix-$1-suffix")
        
    Returns:
        Python-style replacement string (e.g., "\\1" or "prefix-\\1-suffix")
    """
    # First, handle literal $$ (which means a single $ in JavaScript)
    # Replace $$ with a placeholder, then convert backreferences, then restore literal $
    placeholder = '\x00PLACEHOLDER_DOLLAR\x00'
    result = replacement.replace('$$', placeholder)
    
    # Convert $1, $2, $3, etc. to \1, \2, \3, etc.
    # Match $ followed by one or more digits
    def replace_backref(match):
        num = match.group(1)
        return f'\\{num}'  # Use \N format for backreferences
    
    # Replace $N patterns with \N
    result = re.sub(r'\$(\d+)', replace_backref, result)
    
    # Escape any remaining literal $ signs (not part of $N) to $$
    # In Python, $$ in replacement string means a literal $
    result = result.replace('$', '$$')
    
    # Restore placeholder back to $$ (which Python interprets as a single $)
    result = result.replace(placeholder, '$$')
    
    return result


def apply_regex_transform(text: str, regex_pattern: str) -> str:
    """
    Apply a regex transformation to text.
    
    Supports two modes:
    1. Capture group: (.*) - extracts first capture group
    2. Replacement: s/pattern/replacement/ - replaces pattern with replacement
    
    Args:
        text: Text to transform
        regex_pattern: Regex pattern to apply
        
    Returns:
        Transformed text
    """
    if not regex_pattern or regex_pattern == '(.*)':
        return text
    
    # Check for replacement pattern: s/pattern/replacement/
    replacement_match = re.match(r'^s/(.+)/(.*)/$', regex_pattern)
    if replacement_match:
        pattern = replacement_match.group(1)
        replacement = replacement_match.group(2)
        try:
            # Convert JavaScript-style replacement ($1, $2) to Python-style (\1, \2)
            python_replacement = convert_replacement_string(replacement)
            # Use replace without global flag (single replacement)
            result = re.sub(pattern, python_replacement, text, count=1)
            return result
        except re.error as e:
            print(f"Warning: Invalid regex pattern '{regex_pattern}': {e}", file=sys.stderr)
            return text
    
    # Check for capture group pattern
    try:
        match = re.match(regex_pattern, text)
        if match and match.groups():
            return match.group(1)  # Return first capture group
        return text
    except re.error as e:
        print(f"Warning: Invalid regex pattern '{regex_pattern}': {e}", file=sys.stderr)
        return text


def collapse_delimiters(text: str, field_delim: str, subfield_delim: str) -> str:
    """
    Collapse multiple consecutive delimiters and remove leading/trailing delimiters.
    
    Args:
        text: Text to process
        field_delim: Field delimiter
        subfield_delim: Subfield delimiter
        
    Returns:
        Text with collapsed delimiters
    """
    if not text:
        return text
    
    # Escape delimiters for regex
    field_delim_escaped = re.escape(field_delim)
    subfield_delim_escaped = re.escape(subfield_delim)
    
    # Collapse multiple consecutive field delimiters
    result = re.sub(f'{field_delim_escaped}+', field_delim, text)
    
    # Collapse multiple consecutive subfield delimiters within each field
    fields = result.split(field_delim)
    collapsed_fields = []
    for field in fields:
        if field:
            collapsed_field = re.sub(f'{subfield_delim_escaped}+', subfield_delim, field)
            collapsed_fields.append(collapsed_field)
    
    result = field_delim.join(collapsed_fields)
    
    # Remove leading and trailing field delimiters
    leading_delim_regex = re.compile(f'^{field_delim_escaped}+')
    trailing_delim_regex = re.compile(f'{field_delim_escaped}+$')
    result = leading_delim_regex.sub('', result)
    result = trailing_delim_regex.sub('', result)
    
    return result


def transform_filename(
    filename: str,
    filename_regex_map: Dict[int, str],
    current_field_delim: str,
    current_subfield_delim: str,
    new_field_delim: str,
    new_subfield_delim: str
) -> str:
    """
    Apply regex transformations to a filename.
    
    Args:
        filename: Original filename (without extension)
        filename_regex_map: Map of field index to regex pattern
        current_field_delim: Current field delimiter
        current_subfield_delim: Current subfield delimiter
        new_field_delim: New field delimiter for output
        new_subfield_delim: New subfield delimiter for output
        
    Returns:
        Transformed filename
    """
    # Strip file extension
    name_without_ext = re.sub(r'\.(fasta|fa|fas)$', '', filename, flags=re.IGNORECASE)
    
    # Parse into parts (fields and subfields)
    parts = parse_fields(name_without_ext, current_field_delim, current_subfield_delim)
    
    if not parts:
        return name_without_ext
    
    # Apply transformations to each part
    transformed_parts = []
    for part in parts:
        if part['type'] == 'delimiter':
            # Field delimiter - use new field delimiter
            transformed_parts.append(new_field_delim)
        elif part['type'] == 'subdelimiter':
            # Subfield delimiter - use new subfield delimiter
            transformed_parts.append(new_subfield_delim)
        else:
            # Field or subfield - apply transformation
            field_idx = None
            # Find which field this part belongs to by counting field delimiters before it
            field_count = 0
            for p in parts:
                if p == part:
                    break
                if p['type'] == 'delimiter':
                    field_count += 1
            
            # Get regex for this field index
            regex_pattern = filename_regex_map.get(field_count, '(.*)')
            
            # Apply transformation
            transformed = apply_regex_transform(part['text'], regex_pattern)
            
            if transformed:
                transformed_parts.append(transformed)
    
    # Join all parts
    result = ''.join(transformed_parts)
    
    # Amalgamate repeated optional fields
    result = amalgamate_optional_fields(result, new_field_delim, new_subfield_delim)
    
    return result


def amalgamate_optional_fields(text: str, field_delim: str, subfield_delim: str) -> str:
    """
    Amalgamate repeated optional fields in a transformed string.
    Example: "_e-05_e-07" becomes "_e-05-07"
    
    Args:
        text: The transformed text
        field_delim: Field delimiter
        subfield_delim: Subfield delimiter
        
    Returns:
        Text with repeated optional fields amalgamated
    """
    if not text:
        return text
    
    valid_flags = ['a', 'f', 'm', 'w', 'o', 's', 'p', 'e', 'r']
    fields = text.split(field_delim)
    grouped_fields = {}  # flag -> list of values
    result_fields = []
    
    for field in fields:
        if not field:
            # Empty field, keep as is
            result_fields.append(field)
            continue
        
        # Check if this is an optional field (starts with flag-)
        is_optional = False
        flag = None
        value = None
        
        for f in valid_flags:
            if field.startswith(f + subfield_delim):
                is_optional = True
                flag = f
                value = field[len(f) + len(subfield_delim):]
                break
        
        if is_optional:
            # Group this optional field
            if flag not in grouped_fields:
                grouped_fields[flag] = []
            grouped_fields[flag].append(value)
        else:
            # Non-optional field - flush any grouped optional fields before this one
            if grouped_fields:
                for f, values in grouped_fields.items():
                    amalgamated = f + subfield_delim + subfield_delim.join(values)
                    result_fields.append(amalgamated)
                grouped_fields.clear()
            result_fields.append(field)
    
    # Flush any remaining grouped optional fields at the end
    if grouped_fields:
        for f, values in grouped_fields.items():
            amalgamated = f + subfield_delim + subfield_delim.join(values)
            result_fields.append(amalgamated)
    
    return field_delim.join(result_fields)


def transform_sequence_name(
    seq_name: str,
    sequence_regex_map: Dict[int, str],
    current_field_delim: str,
    current_subfield_delim: str,
    new_field_delim: str,
    new_subfield_delim: str
) -> str:
    """
    Apply regex transformations to a sequence name.
    
    Args:
        seq_name: Original sequence name
        sequence_regex_map: Map of field index to regex pattern
        current_field_delim: Current field delimiter
        current_subfield_delim: Current subfield delimiter
        new_field_delim: New field delimiter for output
        new_subfield_delim: New subfield delimiter for output
        
    Returns:
        Transformed sequence name
    """
    # Parse into parts (fields and subfields)
    parts = parse_fields(seq_name, current_field_delim, current_subfield_delim)
    
    if not parts:
        return seq_name
    
    # Apply transformations to each part
    transformed_parts = []
    for part in parts:
        if part['type'] == 'delimiter':
            # Field delimiter - use new field delimiter
            transformed_parts.append(new_field_delim)
        elif part['type'] == 'subdelimiter':
            # Subfield delimiter - use new subfield delimiter
            transformed_parts.append(new_subfield_delim)
        else:
            # Field or subfield - apply transformation
            # Find which field this part belongs to by counting field delimiters before it
            field_count = 0
            for p in parts:
                if p == part:
                    break
                if p['type'] == 'delimiter':
                    field_count += 1
            
            # Get regex for this field index
            regex_pattern = sequence_regex_map.get(field_count, '(.*)')
            
            # Apply transformation
            transformed = apply_regex_transform(part['text'], regex_pattern)
            
            if transformed:
                transformed_parts.append(transformed)
    
    # Join all parts
    result = ''.join(transformed_parts)
    
    # Amalgamate repeated optional fields
    result = amalgamate_optional_fields(result, new_field_delim, new_subfield_delim)
    
    return result


def detect_sequence_type(sequences: List[Tuple[str, str]], skip_count: int = 0) -> str:
    """
    Detect if sequences are amino acid (AA) or nucleotide (NT).
    
    Args:
        sequences: List of (name, sequence) tuples
        skip_count: Number of header sequences to skip
        
    Returns:
        'AA' or 'NT'
    """
    if not sequences or len(sequences) <= skip_count:
        return 'NT'
    
    # Get sample from first non-skipped sequence (up to 100 characters, excluding gaps)
    sample_seq = sequences[skip_count][1]
    sample = ''.join(c for c in sample_seq[:100] if c.upper() not in ['-', 'N', 'X']).upper()
    
    if not sample:
        return 'NT'
    
    # Check for AA-only characters: E, F, I, L, P, Q
    aa_only_chars = set(['E', 'F', 'I', 'L', 'P', 'Q'])
    has_aa_only = any(c in aa_only_chars for c in sample)
    
    if has_aa_only:
        return 'AA'
    
    # Default to NT if no AA-only characters found
    return 'NT'


def validate_transformed_file(
    transformed_filename: str,
    transformed_sequences: List[Tuple[str, str]],
    new_field_delim: str,
    new_subfield_delim: str,
    skip_count: int = 0
) -> Tuple[bool, List[str]]:
    """
    Validate the transformed filename and sequence names.
    
    Args:
        transformed_filename: Transformed filename (without extension)
        transformed_sequences: List of (name, sequence) tuples
        new_field_delim: New field delimiter for collapsing
        new_subfield_delim: New subfield delimiter for collapsing
        skip_count: Number of header sequences to skip from validation
        
    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    errors = []
    
    # Collapse delimiters before validation
    collapsed_filename = collapse_delimiters(transformed_filename, new_field_delim, new_subfield_delim)
    
    # Validate filename
    filename_obj = None
    try:
        filename_obj = parse_fasta_filename(collapsed_filename)
        
        # Check that molecule field matches detected file mode
        if transformed_sequences and len(transformed_sequences) > skip_count:
            detected_mode = detect_sequence_type(transformed_sequences, skip_count).lower()
            filename_mode = filename_obj.molecule.lower()
            if filename_mode != detected_mode:
                errors.append(f"Filename error: Molecule field '{filename_obj.molecule}' does not match detected file type '{detected_mode.upper()}'")
        
        # Check that if filename has a- (alignment) field, sequences are actually aligned
        if filename_obj.alignment and transformed_sequences:
            # Check if all sequences have the same length
            if len(transformed_sequences) > 0:
                first_length = len(transformed_sequences[0][1])
                all_same_length = all(len(seq[1]) == first_length for seq in transformed_sequences)
                
                if not all_same_length:
                    # alignment is a dict with 'main', 'sub', 'subsub' keys
                    alignment_value = filename_obj.alignment.get('main', '') if isinstance(filename_obj.alignment, dict) else str(filename_obj.alignment)
                    errors.append(f"Filename error: Filename contains alignment field 'a-{alignment_value}', but sequences are not aligned (have different lengths)")
    except ValueError as e:
        errors.append(f"Filename error: {e}")
        # Continue to validate sequences even if filename fails
    
    # Validate sequences (skip header sequences)
    for i, (seq_name, seq_data) in enumerate(transformed_sequences[skip_count:], start=skip_count + 1):
        # Collapse delimiters before validation
        collapsed_seq_name = collapse_delimiters(seq_name, new_field_delim, new_subfield_delim)
        
        try:
            seq_name_obj = parse_sequence_name(collapsed_seq_name)
            
            # Check consistency with filename if filename was valid
            if not errors or (errors and 'Filename error' not in errors[0]):
                if filename_obj:
                    if seq_name_obj.network != filename_obj.network:
                        errors.append(f"Sequence {i+1} ('{collapsed_seq_name}'): network mismatch (file: {filename_obj.network}, seq: {seq_name_obj.network})")
                    
                    if seq_name_obj.ptid != filename_obj.ptid:
                        errors.append(f"Sequence {i+1} ('{collapsed_seq_name}'): PTID mismatch (file: {filename_obj.ptid}, seq: {seq_name_obj.ptid})")
        except ValueError as e:
            errors.append(f"Sequence {i+1} ('{collapsed_seq_name}'): {e}")
    
    return len(errors) == 0, errors


def load_json_settings(json_path: Path) -> Dict:
    """
    Load transformation settings from JSON file.
    
    Args:
        json_path: Path to JSON file
        
    Returns:
        Dictionary with settings
    """
    with open(json_path, 'r') as f:
        return json.load(f)


def convert_fasta(
    json_path: Path,
    input_fasta_path: Path,
    output_fasta_path: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    validate: bool = True,
    quiet: bool = False
) -> Tuple[bool, List[str], Optional[str]]:
    """
    Convert a FASTA file to standard compliance using JSON transformation settings.
    
    Args:
        json_path: Path to JSON transformation file
        input_fasta_path: Path to input FASTA file
        output_fasta_path: Path to output FASTA file (default: input name + '_compliant.fasta')
        output_dir: Optional output directory (used if output_fasta_path is None)
        validate: Whether to validate the output
        quiet: Suppress output messages
        
    Returns:
        Tuple of (success, list_of_errors_or_warnings, transformed_filename)
        transformed_filename is the transformed filename (without extension) or None on error
    """
    # Load JSON settings
    try:
        settings = load_json_settings(json_path)
    except Exception as e:
        return False, [f"Error loading JSON file: {e}"], None
    
    # Extract settings
    current_field_delim = settings.get('delimiters', {}).get('current', {}).get('field', '_')
    current_subfield_delim = settings.get('delimiters', {}).get('current', {}).get('subfield', '-')
    new_field_delim = settings.get('delimiters', {}).get('new', {}).get('field', '_')
    new_subfield_delim = settings.get('delimiters', {}).get('new', {}).get('subfield', '-')
    skip_sequences = settings.get('skip_sequences', 0)
    
    # Convert regex maps from string keys to int keys
    filename_transformations = settings.get('transformations', {}).get('filename', {})
    sequence_transformations = settings.get('transformations', {}).get('sequences', {})
    
    filename_regex_map = {int(k): v for k, v in filename_transformations.items()}
    sequence_regex_map = {int(k): v for k, v in sequence_transformations.items()}
    
    # Get input filename (without extension)
    input_filename = input_fasta_path.stem
    
    # Transform filename
    transformed_filename = transform_filename(
        input_filename,
        filename_regex_map,
        current_field_delim,
        current_subfield_delim,
        new_field_delim,
        new_subfield_delim
    )
    
    # Collapse delimiters for output
    collapsed_filename = collapse_delimiters(transformed_filename, new_field_delim, new_subfield_delim)
    
    # Output transformed filename unless quiet
    if not quiet:
        print(f"Transformed filename: {collapsed_filename}")
    
    # Read and transform sequences
    transformed_records = []
    all_transformed_sequences = []
    
    try:
        for i, record in enumerate(SeqIO.parse(input_fasta_path, "fasta")):
            # Use record.description to get the full sequence name (record.id only gets up to first space)
            # Replace spaces in sequence name with current field delimiter (applies to all sequences)
            seq_name = record.description if record.description else record.id
            if seq_name:
                seq_name = re.sub(r'\s+', current_field_delim, seq_name)
            
            # Skip regex transformations for header sequences (first skip_sequences sequences)
            if i < skip_sequences:
                # Keep sequence name with spaces replaced, but no regex transformations
                new_record = SeqRecord(
                    seq=record.seq,
                    id=seq_name,
                    description=""
                )
                transformed_records.append(new_record)
                all_transformed_sequences.append((seq_name, str(record.seq)))
            else:
                # Transform sequence name with regex transformations
                transformed_name = transform_sequence_name(
                    seq_name,
                    sequence_regex_map,
                    current_field_delim,
                    current_subfield_delim,
                    new_field_delim,
                    new_subfield_delim
                )
                
                # Collapse delimiters
                collapsed_name = collapse_delimiters(transformed_name, new_field_delim, new_subfield_delim)
                
                # Create new record with transformed name
                new_record = SeqRecord(
                    seq=record.seq,
                    id=collapsed_name,
                    description=""
                )
                transformed_records.append(new_record)
                all_transformed_sequences.append((collapsed_name, str(record.seq)))
    except Exception as e:
        return False, [f"Error reading FASTA file: {e}"], None
    
    # Validate if requested
    errors = []
    visit_replacement_message = None
    
    if validate:
        is_valid, errors = validate_transformed_file(
            collapsed_filename,
            all_transformed_sequences,
            new_field_delim,
            new_subfield_delim,
            skip_sequences
        )
        
        # If validation succeeds and visit field is "0000", replace it with sorted union of sequence visits
        if is_valid:
            # Parse the transformed filename to get the visit field (field index 2)
            filename_parts = collapsed_filename.split(new_field_delim)
            if len(filename_parts) >= 3:
                visit_field = filename_parts[2]
                
                # Check if visit is "0000"
                if visit_field == '0000':
                    # Collect all visit fields from sequence names (excluding skipped sequences)
                    visit_set = set()
                    sequences_to_check = all_transformed_sequences[skip_sequences:]
                    
                    for seq_name, seq_data in sequences_to_check:
                        seq_parts = seq_name.split(new_field_delim)
                        if len(seq_parts) >= 3:
                            seq_visit = seq_parts[2]
                            # Visit can be a single value or dash-separated (e.g., "0050-0070")
                            # Split by subfield delimiter and add each part
                            visit_parts = seq_visit.split(new_subfield_delim)
                            for v in visit_parts:
                                if v and v.strip():
                                    visit_set.add(v.strip())
                    
                    # Sort visits and join with subfield delimiter
                    sorted_visits = new_subfield_delim.join(sorted(visit_set))
                    
                    if sorted_visits and sorted_visits != '0000':
                        # Replace visit in transformed filename
                        old_visit = filename_parts[2]
                        filename_parts[2] = sorted_visits
                        collapsed_filename = new_field_delim.join(filename_parts)
                        
                        visit_replacement_message = f'Visit field "{old_visit}" was replaced with "{sorted_visits}" (sorted union of visits from sequence names).'
                        
                        # Re-validate with updated filename
                        is_valid, errors = validate_transformed_file(
                            collapsed_filename,
                            all_transformed_sequences,
                            new_field_delim,
                            new_subfield_delim,
                            skip_sequences
                        )
        
        if not is_valid:
            print("Validation errors found:", file=sys.stderr)
            for error in errors:
                print(f"  - {error}", file=sys.stderr)
            if visit_replacement_message:
                print(f"Note: {visit_replacement_message}", file=sys.stderr)
            return False, errors
        
        if visit_replacement_message:
            print(f"Note: {visit_replacement_message}")
    
    # Determine output path
    if output_fasta_path is None:
        if output_dir is not None:
            # Use output directory with transformed filename
            output_fasta_path = output_dir / f"{collapsed_filename}_compliant.fasta"
        else:
            # Use transformed filename with _compliant suffix in input directory
            output_fasta_path = input_fasta_path.parent / f"{collapsed_filename}_compliant.fasta"
    
    # Write output FASTA file with 80 characters per line (matching compliance app)
    try:
        with open(output_fasta_path, 'w') as f:
            for record in transformed_records:
                # Write header line
                f.write(f">{record.id}\n")
                # Write sequence in 80-character lines
                seq_str = str(record.seq)
                for i in range(0, len(seq_str), 80):
                    f.write(seq_str[i:i+80] + "\n")
        print(f"✓ Converted FASTA file written to: {output_fasta_path}")
        if validate and not errors:
            print("✓ Validation passed: File is compliant with standard convention")
    except Exception as e:
        return False, [f"Error writing output file: {e}"], None
    
    return True, errors if errors else ["Conversion completed successfully"], collapsed_filename


def main():
    """Main entry point."""
    if len(sys.argv) < 3:
        print("Usage (single file): python fasta2standard.py <json_file> <input_fasta> [output_fasta] [-quiet]", file=sys.stderr)
        print("Usage (directory): python fasta2standard.py <json_file> <input_directory> [-quiet]", file=sys.stderr)
        sys.exit(1)
    
    # Parse arguments
    args = sys.argv[1:]
    quiet = '-quiet' in args
    if quiet:
        args.remove('-quiet')
    
    json_path = Path(args[0])
    input_path = Path(args[1])
    output_path = Path(args[2]) if len(args) > 2 else None
    
    # Validate JSON file exists
    if not json_path.exists():
        print(f"Error: JSON file not found: {json_path}", file=sys.stderr)
        sys.exit(1)
    
    if not input_path.exists():
        print(f"Error: Input path not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    # Check if input is a directory
    if input_path.is_dir():
        # Process directory
        if not quiet:
            print(f"Processing directory: {input_path}")
        
        # Create output directory with _compliant suffix
        output_dir = Path(str(input_path) + "_compliant")
        output_dir.mkdir(exist_ok=True)
        
        if not quiet:
            print(f"Output directory: {output_dir}")
        
        # Find all FASTA files in the directory
        fasta_extensions = ['*.fasta', '*.fa', '*.fas', '*.FASTA', '*.FA', '*.FAS']
        fasta_files = []
        for ext in fasta_extensions:
            fasta_files.extend(input_path.glob(ext))
        
        if not fasta_files:
            print(f"Warning: No FASTA files found in {input_path}", file=sys.stderr)
            sys.exit(1)
        
        # Process each FASTA file
        success_count = 0
        error_count = 0
        
        for fasta_file in sorted(fasta_files):
            if not quiet:
                print(f"\nProcessing: {fasta_file.name}")
            
            # Call convert_fasta with output directory to use transformed filename
            success, messages, transformed_filename = convert_fasta(
                json_path, 
                fasta_file, 
                None,  # Let convert_fasta determine the transformed filename
                output_dir=output_dir,  # Specify output directory
                validate=True, 
                quiet=quiet
            )
            
            if success and transformed_filename and not quiet:
                print(f"  Output: {transformed_filename}_compliant.fasta")
            
            if success:
                success_count += 1
                if not quiet:
                    for msg in messages:
                        if "error" in msg.lower() or "warning" in msg.lower():
                            print(f"  {msg}")
            else:
                error_count += 1
                print(f"Error processing {fasta_file.name}:", file=sys.stderr)
                for msg in messages:
                    print(f"  {msg}", file=sys.stderr)
        
        # Summary
        if not quiet:
            print(f"\nSummary: {success_count} files processed successfully, {error_count} errors")
        
        if error_count > 0:
            sys.exit(1)
    else:
        # Process single file
        # Convert
        success, messages, _ = convert_fasta(json_path, input_path, output_path, validate=True, quiet=quiet)
        
        if not success:
            sys.exit(1)


if __name__ == "__main__":
    main()

