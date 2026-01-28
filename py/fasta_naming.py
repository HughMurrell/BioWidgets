#!/usr/bin/env python3
"""
FASTA File and Sequence Naming Convention Utility

This standalone script provides parsing, validation, and conversion utilities
for FASTA files following the network/protocol naming convention.

File naming format:
<network><protocol>_<PTID>_<visit>_<region>_<molecule>_[optional fields].fasta

Sequence naming format:
<network><cohort>_<PTID>_<visit>_<region>_<sequencing method>-<sequence ID>_[optional fields]

Usage:
    # As a module
    from fasta_naming import parse_fasta_filename, FastaFileNameBuilder
    
    # From command line
    python fasta_naming.py validate myfile.fasta
    python fasta_naming.py parse myfile.fasta

Author: Chris Barry
Date: 2025-11-24
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional, Union
import re
import sys

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython", file=sys.stderr)
    sys.exit(1)


# ============================================================================
# Global Validation Constants
# ============================================================================

# Regex patterns for validation
NETWORK_PATTERN = r'^[A-Z]$'
PROTOCOL_PATTERN = r'^\d{3,4}[A-Z]?$'
PTID_PATTERN = r'^(CAP\d{3}|\d{4,6})$'
VISIT_PATTERN = r'^\d{4}(-\d{4})*$'
SEQUENCE_ID_PATTERN = r'^\w{8}$'

# Valid values
VALID_REGIONS = ['env', 'pol', 'gag', 'nef', 'ren', 'gp']
VALID_MOLECULES = ['aa', 'nt']
VALID_SEQUENCING_METHODS = ['pb', 'sa', 'np', 'il', 'sy']

# Optional field flags and their corresponding attribute names
OPTIONAL_FIELD_FLAGS = {
    'a': 'alignment',
    'f': 'filters',
    'm': 'modifiers',
    'w': 'additions',
    'o': 'removals',
    's': 'sequencing',
    'p': 'processing',
    'e': 'extras',
    'r': 'reviewed'
}


# ============================================================================
# Data Classes for Structured Naming
# ============================================================================

def _format_optional_field(flag: str, field_data: Dict[str, any]) -> str:
    """
    Format optional field dictionary into string representation.
    
    Args:
        flag: Single-letter flag (a, f, m, etc.)
        field_data: Dict with 'main', 'sub', 'subsub' keys
        
    Returns:
        Formatted string like "a-maff" or "a-pwis-HXB2+ref" or "m-rw+CB+AG"
        
    Examples:
        {'main': 'maff', 'sub': [], 'subsub': {}} -> "a-maff"
        {'main': 'pwis', 'sub': ['HXB2'], 'subsub': {}} -> "a-pwis-HXB2"
        {'main': 'rw', 'sub': [], 'subsub': {'rw': ['CB', 'AG']}} -> "m-rw+CB+AG"
        {'main': 'method', 'sub': ['s1', 's2'], 'subsub': {'s1': ['x', 'y']}} -> "a-method-s1+x+y-s2"
    """
    parts = [f"{flag}-{field_data['main']}"]
    
    # Add subcategories
    for sub in field_data.get('sub', []):
        if sub in field_data.get('subsub', {}):
            # Has subsubcategories
            subsub_parts = '+'.join(field_data['subsub'][sub])
            parts.append(f"{sub}+{subsub_parts}")
        else:
            # No subsubcategories
            parts.append(sub)
    
    # Add subsub entries that don't have explicit sub entries
    for key, values in field_data.get('subsub', {}).items():
        if key == field_data['main'] and key not in field_data.get('sub', []):
            # subsubcategories directly attached to main
            subsub_str = '+'.join(values)
            parts[0] += f"+{subsub_str}"
    
    return '-'.join(parts)


def _parse_optional_field(field_str: str) -> Dict[str, any]:
    """
    Parse optional field string into structured dictionary.
    
    Args:
        field_str: String like "a-maff" or "f-func-ctm1" or "m-rw+CB+AG"
        
    Returns:
        Dict with 'main', 'sub', 'subsub' keys
        
    Examples:
        "a-maff" -> {'main': 'maff', 'sub': [], 'subsub': {}}
        "a-pwis-HXB2" -> {'main': 'pwis', 'sub': ['HXB2'], 'subsub': {}}
        "m-rw+CB+AG" -> {'main': 'rw', 'sub': [], 'subsub': {'rw': ['CB', 'AG']}}
        "a-method-s1+x+y-s2" -> {'main': 'method', 'sub': ['s1', 's2'], 'subsub': {'s1': ['x', 'y']}}
    """
    # Remove flag prefix (e.g., "a-" or "f-")
    if '-' not in field_str:
        raise ValueError(f"Invalid optional field format: {field_str}")
    
    content = field_str.split('-', 1)[1]  # Remove flag, keep everything after first dash
    
    # Split by dash to get main and subcategories
    parts = content.split('-')
    
    result = {
        'main': '',
        'sub': [],
        'subsub': {}
    }
    
    # First part might contain main+subsub (e.g., "rw+CB+AG")
    if '+' in parts[0]:
        main_parts = parts[0].split('+')
        result['main'] = main_parts[0]
        result['subsub'][main_parts[0]] = main_parts[1:]
    else:
        result['main'] = parts[0]
    
    # Process remaining parts (subcategories)
    for part in parts[1:]:
        if '+' in part:
            # Subcategory with subsubcategories
            sub_parts = part.split('+')
            sub_name = sub_parts[0]
            result['sub'].append(sub_name)
            result['subsub'][sub_name] = sub_parts[1:]
        else:
            # Simple subcategory
            result['sub'].append(part)
    
    return result


@dataclass
class FastaFileName:
    """Parsed FASTA filename following naming convention."""
    
    # Obligate fields
    network: str  # 1 letter: V, C, U, Z
    protocol: str  # 3-4 chars: 804, 705, 012C
    ptid: str  # 4-6 chars: 0123, CAP001
    visit: str  # 4 digits: 0000, 1000, 2000-3000
    region: str  # 3 letters: env, pol, gag, nef, gp
    molecule: str  # 2 letters: aa, nt
    
    # Optional fields stored as structured data
    # Format: {'main': str, 'sub': list, 'subsub': dict}
    # Example: {'main': 'maff', 'sub': ['method1', 'method2'], 'subsub': {'method1': ['HXB2', 'ref']}}
    alignment: Optional[Dict[str, any]] = None  # a-maff, a-pwis+HXB2
    filters: Optional[Dict[str, any]] = None  # f-func, f-ctm1, f-drep
    modifiers: Optional[Dict[str, any]] = None  # m-cent, m-rw+CB
    additions: Optional[Dict[str, any]] = None  # w-HXB2, w-ref
    removals: Optional[Dict[str, any]] = None  # o-seq123
    sequencing: Optional[Dict[str, any]] = None  # s-2024-p01-nicd
    processing: Optional[Dict[str, any]] = None  # p-v1p5-qt97p5
    extras: Optional[Dict[str, any]] = None  # e-custom
    reviewed: Optional[Dict[str, any]] = None  # r-CB, r-AG
    
    def __post_init__(self):
        """Validate fields after initialization."""
        self.validate()
    
    def validate(self) -> None:
        """Validate all fields according to convention rules."""
        # Network validation
        if not re.match(NETWORK_PATTERN, self.network):
            raise ValueError(f"Invalid network '{self.network}': must be 1 uppercase letter")
        
        # Protocol validation
        if not re.match(PROTOCOL_PATTERN, self.protocol):
            raise ValueError(f"Invalid protocol '{self.protocol}': must be 3-4 characters")
        
        # PTID validation
        if not re.match(PTID_PATTERN, self.ptid):
            raise ValueError(f"Invalid PTID '{self.ptid}': must be 4-6 digits or CAP### format")
        
        # Visit validation
        if not re.match(VISIT_PATTERN, self.visit):
            raise ValueError(f"Invalid visit '{self.visit}': must be 4 digits or dash-separated")
        
        # Region validation
        if self.region.lower() not in VALID_REGIONS:
            raise ValueError(f"Invalid region '{self.region}': must be one of {VALID_REGIONS}")
        
        # Molecule validation
        if self.molecule.lower() not in VALID_MOLECULES:
            raise ValueError(f"Invalid molecule '{self.molecule}': must be one of {VALID_MOLECULES}")
    
    def to_filename(self) -> str:
        """Generate filename string from components."""
        parts = [
            f"{self.network}{self.protocol}",
            self.ptid,
            self.visit,
            self.region,
            self.molecule
        ]
        
        # Add optional fields if present
        if self.alignment:
            parts.append(_format_optional_field('a', self.alignment))
        if self.filters:
            parts.append(_format_optional_field('f', self.filters))
        if self.modifiers:
            parts.append(_format_optional_field('m', self.modifiers))
        if self.additions:
            parts.append(_format_optional_field('w', self.additions))
        if self.removals:
            parts.append(_format_optional_field('o', self.removals))
        if self.sequencing:
            parts.append(_format_optional_field('s', self.sequencing))
        if self.processing:
            parts.append(_format_optional_field('p', self.processing))
        if self.extras:
            parts.append(_format_optional_field('e', self.extras))
        if self.reviewed:
            parts.append(_format_optional_field('r', self.reviewed))
        
        return "_".join(parts) + ".fasta"
    
    def get_optional_field(self, field_name: str, subcategory: Optional[str] = None) -> any:
        """
        Query optional field data.
        
        Args:
            field_name: Name of optional field ('alignment', 'filters', etc.)
            subcategory: Optional subcategory name to retrieve specific data
            
        Returns:
            Full field dict if no subcategory, specific subcategory data otherwise
            
        Examples:
            fn.get_optional_field('alignment')  # Returns full alignment dict
            fn.get_optional_field('alignment', 'main')  # Returns 'maff'
            fn.get_optional_field('alignment', 'subsub')  # Returns subsub dict
        """
        field_data = getattr(self, field_name, None)
        if field_data is None:
            return None
        if subcategory:
            return field_data.get(subcategory)
        return field_data
    
    def set_optional_field(self, field_name: str, main: str, 
                          sub: Optional[List[str]] = None,
                          subsub: Optional[Dict[str, List[str]]] = None) -> None:
        """
        Set optional field with structured data.
        
        Args:
            field_name: Name of optional field ('alignment', 'filters', etc.)
            main: Main category value
            sub: List of subcategories (optional)
            subsub: Dict mapping subcategories to subsubcategories (optional)
            
        Examples:
            fn.set_optional_field('alignment', 'maff')
            fn.set_optional_field('alignment', 'pwis', sub=['HXB2'])
            fn.set_optional_field('modifiers', 'rw', subsub={'rw': ['CB', 'AG']})
        """
        setattr(self, field_name, {
            'main': main,
            'sub': sub or [],
            'subsub': subsub or {}
        })
    
    def add_subcategory(self, field_name: str, subcategory: str,
                       subsubcategories: Optional[List[str]] = None) -> None:
        """
        Add a subcategory to an existing optional field.
        
        Args:
            field_name: Name of optional field
            subcategory: Subcategory to add
            subsubcategories: Optional list of subsubcategories for this subcategory
            
        Example:
            fn.add_subcategory('alignment', 'HXB2')
            fn.add_subcategory('modifiers', 'rw', ['CB', 'AG'])
        """
        field_data = getattr(self, field_name, None)
        if field_data is None:
            raise ValueError(f"Field '{field_name}' is not set. Use set_optional_field first.")
        
        if subcategory not in field_data['sub']:
            field_data['sub'].append(subcategory)
        
        if subsubcategories:
            field_data['subsub'][subcategory] = subsubcategories


@dataclass
class SequenceName:
    """Parsed sequence name following naming convention."""
    
    # Obligate fields
    network: str  # 1 letter
    cohort: str  # 3-4 chars
    ptid: str  # 6 chars
    visit: str  # 4 digits
    region: str  # 3 letters
    sequencing_method: str  # 2 letters: pb, sa, np
    sequence_id: str  # 8 digits (or empty if collapsed)
    
    # Optional fields
    processing: Optional[str] = None  # p-fs45-ma90
    modifiers: Optional[str] = None  # m-mn+CB
    extras: Optional[str] = None  # e-custom
    collapsed: Optional[str] = None  # coll-1-523
    
    def __post_init__(self):
        """Validate fields after initialization."""
        self.validate()
    
    def validate(self) -> None:
        """Validate sequence name fields."""
        # Network validation
        if not re.match(NETWORK_PATTERN, self.network):
            raise ValueError(f"Invalid network '{self.network}': must be 1 uppercase letter")
        
        # Sequencing method validation
        if self.sequencing_method.lower() not in VALID_SEQUENCING_METHODS:
            raise ValueError(f"Invalid sequencing method '{self.sequencing_method}': must be one of {VALID_SEQUENCING_METHODS}")
        
        # Sequence ID validation (optional for collapsed sequences)
        if self.sequence_id and not self.collapsed:
            if not re.match(SEQUENCE_ID_PATTERN, self.sequence_id):
                raise ValueError(f"Invalid sequence ID '{self.sequence_id}': must be 8 characters")
    
    def to_sequence_id(self) -> str:
        """Generate sequence ID string from components."""
        parts = [
            f"{self.network}{self.cohort}",
            self.ptid,
            self.visit,
            self.region,
            f"{self.sequencing_method}-{self.sequence_id}" if self.sequence_id else self.sequencing_method
        ]
        
        # Add optional fields
        if self.collapsed:
            parts.append(self.collapsed)
        if self.processing:
            parts.append(self.processing)
        if self.modifiers:
            parts.append(self.modifiers)
        if self.extras:
            parts.append(self.extras)
        
        return "_".join(parts)


# ============================================================================
# Parser Functions
# ============================================================================

def parse_fasta_filename(filename: str) -> FastaFileName:
    """
    Parse a FASTA filename into structured components.
    
    Args:
        filename: FASTA filename (with or without .fasta extension)
        
    Returns:
        FastaFileName object with parsed components
        
    Raises:
        ValueError: If filename doesn't match convention
    """
    # Remove .fasta extension if present
    name = filename.replace('.fasta', '')
    
    # Split by underscore
    parts = name.split('_')
    
    if len(parts) < 5:
        raise ValueError(f"Filename has too few parts: expected at least 5, got {len(parts)}")
    
    # Parse obligate fields
    network_protocol = parts[0]
    network = network_protocol[0]
    protocol = network_protocol[1:]
    
    ptid = parts[1]
    visit = parts[2]
    region = parts[3]
    molecule = parts[4]
    
    # Parse optional fields by flag
    optional_fields = {}
    for part in parts[5:]:
        if part.startswith('a-'):
            optional_fields['alignment'] = _parse_optional_field(part)
        elif part.startswith('f-'):
            optional_fields['filters'] = _parse_optional_field(part)
        elif part.startswith('m-'):
            optional_fields['modifiers'] = _parse_optional_field(part)
        elif part.startswith('w-'):
            optional_fields['additions'] = _parse_optional_field(part)
        elif part.startswith('o-'):
            optional_fields['removals'] = _parse_optional_field(part)
        elif part.startswith('s-'):
            optional_fields['sequencing'] = _parse_optional_field(part)
        elif part.startswith('p-'):
            optional_fields['processing'] = _parse_optional_field(part)
        elif part.startswith('e-'):
            optional_fields['extras'] = _parse_optional_field(part)
        elif part.startswith('r-'):
            optional_fields['reviewed'] = _parse_optional_field(part)
    
    return FastaFileName(
        network=network,
        protocol=protocol,
        ptid=ptid,
        visit=visit,
        region=region,
        molecule=molecule,
        **optional_fields
    )


def parse_sequence_name(seq_id: str) -> SequenceName:
    """
    Parse a sequence ID into structured components.
    
    Args:
        seq_id: Sequence identifier from FASTA file
        
    Returns:
        SequenceName object with parsed components
        
    Raises:
        ValueError: If sequence ID doesn't match convention
    """
    parts = seq_id.split('_')
    
    if len(parts) < 5:
        raise ValueError(f"Sequence ID has too few parts: expected at least 5, got {len(parts)}")
    
    # Parse obligate fields
    network_cohort = parts[0]
    network = network_cohort[0]
    cohort = network_cohort[1:]
    
    ptid = parts[1]
    visit = parts[2]
    region = parts[3]
    
    # Parse sequencing method and ID
    seq_parts = parts[4].split('-')
    sequencing_method = seq_parts[0]
    sequence_id = seq_parts[1] if len(seq_parts) > 1 else ""
    
    # Parse optional fields
    optional_fields = {}
    for part in parts[5:]:
        if part.startswith('p-'):
            optional_fields['processing'] = part
        elif part.startswith('m-'):
            optional_fields['modifiers'] = part
        elif part.startswith('e-'):
            optional_fields['extras'] = part
        elif part.startswith('coll-'):
            optional_fields['collapsed'] = part
    
    return SequenceName(
        network=network,
        cohort=cohort,
        ptid=ptid,
        visit=visit,
        region=region,
        sequencing_method=sequencing_method,
        sequence_id=sequence_id,
        **optional_fields
    )


# ============================================================================
# Builder Classes for User-Friendly Input
# ============================================================================

class FastaFileNameBuilder:
    """Builder for constructing FASTA filenames with validation."""
    
    def __init__(self):
        self._data = {}
    
    def set_network(self, network: str) -> 'FastaFileNameBuilder':
        """Set network (1 letter: V, C, U, Z)."""
        self._data['network'] = network.upper()
        return self
    
    def set_protocol(self, protocol: str) -> 'FastaFileNameBuilder':
        """Set protocol (3-4 chars: 804, 705, 012C)."""
        self._data['protocol'] = protocol
        return self
    
    def set_ptid(self, ptid: str) -> 'FastaFileNameBuilder':
        """Set participant ID (4-6 chars or CAP###)."""
        self._data['ptid'] = ptid
        return self
    
    def set_visit(self, visit: Union[str, int]) -> 'FastaFileNameBuilder':
        """Set visit (4 digits or dash-separated)."""
        if isinstance(visit, int):
            visit = f"{visit:04d}"
        self._data['visit'] = visit
        return self
    
    def set_region(self, region: str) -> 'FastaFileNameBuilder':
        """Set gene region (env, pol, gag, nef, ren, gp)."""
        self._data['region'] = region.lower()
        return self
    
    def set_molecule(self, molecule: str) -> 'FastaFileNameBuilder':
        """Set molecule type (aa, nt)."""
        self._data['molecule'] = molecule.lower()
        return self
    
    def set_alignment(self, method: str, params: Optional[str] = None) -> 'FastaFileNameBuilder':
        """Set alignment method (e.g., 'maff', 'pwis')."""
        self._data['alignment'] = {
            'main': method,
            'sub': [params] if params else [],
            'subsub': {}
        }
        return self
    
    def set_filter(self, filter_type: str) -> 'FastaFileNameBuilder':
        """Set filter type (e.g., 'func', 'ctm1', 'drep')."""
        self._data['filters'] = {
            'main': filter_type,
            'sub': [],
            'subsub': {}
        }
        return self
    
    def set_modifier(self, modifier: str, *submodifiers: str) -> 'FastaFileNameBuilder':
        """
        Set modifier (e.g., 'cent', 'rw').
        
        Args:
            modifier: Main modifier
            *submodifiers: Optional submodifiers (e.g., 'CB', 'AG')
            
        Example:
            builder.set_modifier('rw', 'CB', 'AG')  # Results in m-rw+CB+AG
        """
        self._data['modifiers'] = {
            'main': modifier,
            'sub': [],
            'subsub': {modifier: list(submodifiers)} if submodifiers else {}
        }
        return self
    
    def add_reference(self, ref: str) -> 'FastaFileNameBuilder':
        """Add reference sequence (e.g., 'HXB2', 'ref')."""
        self._data['additions'] = {
            'main': ref,
            'sub': [],
            'subsub': {}
        }
        return self
    
    def set_removal(self, removal: str) -> 'FastaFileNameBuilder':
        """Set removed sequence identifier."""
        self._data['removals'] = {
            'main': removal,
            'sub': [],
            'subsub': {}
        }
        return self
    
    def set_sequencing(self, year: int, pool: int, facility: str) -> 'FastaFileNameBuilder':
        """Set sequencing details."""
        self._data['sequencing'] = {
            'main': str(year),
            'sub': [f"p{pool:02d}", facility.lower()],
            'subsub': {}
        }
        return self
    
    def set_processing(self, version: Optional[str] = None, **params) -> 'FastaFileNameBuilder':
        """Set processing parameters."""
        parts = []
        if version:
            parts.append(f"v{version}")
        for key, value in params.items():
            parts.append(f"{key}{value}")
        
        self._data['processing'] = {
            'main': parts[0] if parts else '',
            'sub': parts[1:] if len(parts) > 1 else [],
            'subsub': {}
        }
        return self
    
    def set_extras(self, extras: str) -> 'FastaFileNameBuilder':
        """Set extra information."""
        self._data['extras'] = {
            'main': extras,
            'sub': [],
            'subsub': {}
        }
        return self
    
    def set_reviewer(self, *initials: str) -> 'FastaFileNameBuilder':
        """
        Set reviewer initials.
        
        Args:
            *initials: One or more reviewer initials
            
        Example:
            builder.set_reviewer('CB')  # Results in r-CB
            builder.set_reviewer('CB', 'AG')  # Results in r-CB+AG
        """
        if not initials:
            raise ValueError("At least one reviewer initial required")
        
        if len(initials) == 1:
            self._data['reviewed'] = {
                'main': initials[0].upper(),
                'sub': [],
                'subsub': {}
            }
        else:
            self._data['reviewed'] = {
                'main': initials[0].upper(),
                'sub': [],
                'subsub': {initials[0].upper(): [i.upper() for i in initials[1:]]}
            }
        return self
    
    def build(self) -> FastaFileName:
        """Build and validate the FastaFileName object."""
        required = ['network', 'protocol', 'ptid', 'visit', 'region', 'molecule']
        missing = [f for f in required if f not in self._data]
        if missing:
            raise ValueError(f"Missing required fields: {', '.join(missing)}")
        
        return FastaFileName(**self._data)


class SequenceNameBuilder:
    """Builder for constructing sequence names with validation."""
    
    def __init__(self):
        self._data = {}
    
    def set_network(self, network: str) -> 'SequenceNameBuilder':
        """Set network (1 letter)."""
        self._data['network'] = network.upper()
        return self
    
    def set_cohort(self, cohort: str) -> 'SequenceNameBuilder':
        """Set cohort (3-4 chars)."""
        self._data['cohort'] = cohort
        return self
    
    def set_ptid(self, ptid: str) -> 'SequenceNameBuilder':
        """Set participant ID."""
        self._data['ptid'] = ptid
        return self
    
    def set_visit(self, visit: Union[str, int]) -> 'SequenceNameBuilder':
        """Set visit (4 digits)."""
        if isinstance(visit, int):
            visit = f"{visit:04d}"
        self._data['visit'] = visit
        return self
    
    def set_region(self, region: str) -> 'SequenceNameBuilder':
        """Set gene region."""
        self._data['region'] = region.lower()
        return self
    
    def set_sequencing(self, method: str, seq_id: str) -> 'SequenceNameBuilder':
        """Set sequencing method and ID."""
        self._data['sequencing_method'] = method.lower()
        self._data['sequence_id'] = seq_id
        return self
    
    def set_collapsed(self, rank: int, count: int) -> 'SequenceNameBuilder':
        """Set collapsed sequence info."""
        self._data['collapsed'] = f"coll-{rank}-{count}"
        self._data['sequence_id'] = ""  # Collapsed sequences don't have IDs
        return self
    
    def set_processing(self, **params) -> 'SequenceNameBuilder':
        """Set processing parameters (e.g., fs=45, ma=90)."""
        parts = [f"{k}{v}" for k, v in params.items()]
        self._data['processing'] = f"p-{'-'.join(parts)}"
        return self
    
    def set_modifier(self, modifier: str) -> 'SequenceNameBuilder':
        """Set modifier (e.g., 'mn+CB', 'rw+AG')."""
        self._data['modifiers'] = f"m-{modifier}"
        return self
    
    def set_extras(self, extras: str) -> 'SequenceNameBuilder':
        """Set extra information."""
        self._data['extras'] = f"e-{extras}"
        return self
    
    def build(self) -> SequenceName:
        """Build and validate the SequenceName object."""
        required = ['network', 'cohort', 'ptid', 'visit', 'region', 'sequencing_method']
        missing = [f for f in required if f not in self._data]
        if missing:
            raise ValueError(f"Missing required fields: {', '.join(missing)}")
        
        if 'sequence_id' not in self._data:
            self._data['sequence_id'] = ""
        
        return SequenceName(**self._data)


# ============================================================================
# File I/O Functions
# ============================================================================

def read_fasta_with_naming(fasta_path: Path) -> Dict[str, any]:
    """
    Read FASTA file and parse naming conventions.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Dictionary with 'filename', 'sequences' keys
    """
    filename_obj = parse_fasta_filename(fasta_path.name)
    
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        try:
            seq_name = parse_sequence_name(record.id)
            sequences.append({
                'record': record,
                'parsed_name': seq_name
            })
        except ValueError as e:
            print(f"Warning: Could not parse sequence name '{record.id}': {e}", file=sys.stderr)
            sequences.append({
                'record': record,
                'parsed_name': None
            })
    
    return {
        'filename': filename_obj,
        'sequences': sequences
    }


def _parse_filename_pattern(filename: str, pattern: str) -> Dict[str, str]:
    """
    Parse a filename using a pattern template with category placeholders.
    
    Args:
        filename: The actual filename to parse
        pattern: Pattern with categories in <brackets> and literal characters
                Example: "<network><protocol>_<ptid>_<visit>_<region>_<molecule>.fasta"
                Example: "<ptid>-<visit>.fasta"
        
    Returns:
        Dictionary mapping category names to extracted values
        
    Raises:
        ValueError: If filename doesn't match pattern
        
    Note:
        Patterns work best with clear literal separators between categories.
        Avoid patterns with adjacent placeholders without separators or with
        multiple consecutive separator characters (e.g., "----") between placeholders.
    """
    # Split pattern into parts (categories and literals)
    import re
    
    # Find all categories and their positions
    category_pattern = r'<([^>]+)>'
    parts = []
    last_end = 0
    
    for match in re.finditer(category_pattern, pattern):
        # Add literal before this category
        if match.start() > last_end:
            parts.append(('literal', pattern[last_end:match.start()]))
        # Add this category
        parts.append(('category', match.group(1)))
        last_end = match.end()
    
    # Add trailing literal
    if last_end < len(pattern):
        parts.append(('literal', pattern[last_end:]))
    
    # Build regex pattern
    regex_parts = ['^']
    category_mapping = {}
    
    for i, (part_type, value) in enumerate(parts):
        if part_type == 'literal':
            # Escape and add literal
            regex_parts.append(re.escape(value))
        else:
            # Category - sanitize name and determine what follows
            sanitized = value.replace('-', '_')
            category_mapping[sanitized] = value
            
            # Look ahead to see what comes next
            if i + 1 < len(parts) and parts[i + 1][0] == 'literal':
                # Next is a literal - match until that literal
                next_literal = parts[i + 1][1]
                if next_literal:
                    # Use non-greedy match that stops before the next literal
                    regex_parts.append(f"(?P<{sanitized}>.+?)(?={re.escape(next_literal)})")
                else:
                    # Empty literal (shouldn't happen but handle it)
                    regex_parts.append(f"(?P<{sanitized}>.+?)")
            else:
                # No literal follows (end of pattern or another category follows)
                # Use greedy match to consume as much as possible
                regex_parts.append(f"(?P<{sanitized}>.+?)")
    
    regex_parts.append('$')
    regex_pattern = ''.join(regex_parts)
    
    # Match filename against pattern
    match = re.match(regex_pattern, filename)
    if not match:
        raise ValueError(f"Filename '{filename}' doesn't match pattern '{pattern}'")
    
    # Convert sanitized names back to original category names
    result = {}
    for sanitized, value in match.groupdict().items():
        original = category_mapping.get(sanitized, sanitized)
        result[original] = value
    
    return result


def load_fasta_with_metadata(
    fasta_path: Path,
    pattern: Optional[str] = None,
    network: Optional[str] = None,
    protocol: Optional[str] = None,
    ptid: Optional[str] = None,
    visit: Optional[Union[str, int]] = None,
    region: Optional[str] = None,
    molecule: Optional[str] = None,
    **optional_fields
) -> Dict[str, any]:
    """
    Load FASTA file and assign naming convention metadata on-the-fly.
    
    This function allows you to load any FASTA file (even with non-conforming 
    filenames) and assign naming convention metadata to it. Useful for converting
    legacy files or files from external sources.
    
    Args:
        fasta_path: Path to FASTA file
        pattern: Optional filename pattern with categories in <brackets>.
                Example: "<network><protocol>_<ptid>_<visit>_<region>_<molecule>.fasta"
                Example: "<ptid>-<visit>----<region>.fa"
                Categories are extracted and used as metadata fields.
        network: Network identifier (1 letter: V, C, U, Z)
        protocol: Protocol identifier (3-4 chars: 804, 705, 012C)
        ptid: Participant ID (4-6 chars or CAP###)
        visit: Visit identifier (4 digits or dash-separated)
        region: Gene region (env, pol, gag, nef, ren, gp)
        molecule: Molecule type (aa, nt)
        **optional_fields: Optional fields like alignment, filters, modifiers, etc.
                          Use 'alignment_method', 'filter_type', 'modifier', 
                          'reference', 'reviewer', etc.
        
    Returns:
        Dictionary with 'filename', 'sequences', 'original_filename' keys
        
    Examples:
        # Load non-conforming file with metadata
        data = load_fasta_with_metadata(
            Path("my_sequences.fasta"),
            network='V',
            protocol='804',
            ptid='CAP001',
            visit=0,
            region='env',
            molecule='nt'
        )
        
        # Use pattern to parse filename
        data = load_fasta_with_metadata(
            Path("CAP001-0000----env.fasta"),
            pattern="<ptid>-<visit>----<region>.fasta",
            network='V',
            protocol='804',
            molecule='nt'
        )
        
        # Pattern with complex separators
        data = load_fasta_with_metadata(
            Path("V804_CAP001_visit0000_env_nt.fasta"),
            pattern="<network><protocol>_<ptid>_visit<visit>_<region>_<molecule>.fasta"
        )
    """
    # Try to parse existing filename first
    try:
        return read_fasta_with_naming(fasta_path)
    except ValueError:
        pass  # Filename doesn't conform, continue with manual metadata
    
    # If pattern provided, extract metadata from filename
    extracted_fields = {}
    if pattern:
        try:
            extracted_fields = _parse_filename_pattern(fasta_path.name, pattern)
        except ValueError as e:
            raise ValueError(f"Pattern parsing failed: {e}")
    
    # Merge extracted fields with explicitly provided fields (explicit takes precedence)
    merged_fields = {
        'network': network if network is not None else extracted_fields.get('network'),
        'protocol': protocol if protocol is not None else extracted_fields.get('protocol'),
        'ptid': ptid if ptid is not None else extracted_fields.get('ptid'),
        'visit': visit if visit is not None else extracted_fields.get('visit'),
        'region': region if region is not None else extracted_fields.get('region'),
        'molecule': molecule if molecule is not None else extracted_fields.get('molecule')
    }
    
    # Check if all required fields are provided
    required_fields = merged_fields
    
    missing = [k for k, v in required_fields.items() if v is None]
    if missing:
        raise ValueError(
            f"Filename '{fasta_path.name}' doesn't conform to naming convention. "
            f"Please provide required fields: {', '.join(missing)}"
        )
    
    # Build filename object from merged metadata
    builder = FastaFileNameBuilder()
    builder.set_network(merged_fields['network'])
    builder.set_protocol(merged_fields['protocol'])
    builder.set_ptid(merged_fields['ptid'])
    builder.set_visit(merged_fields['visit'])
    builder.set_region(merged_fields['region'])
    builder.set_molecule(merged_fields['molecule'])
    
    # Handle optional fields with semantic parameter names
    if 'alignment_method' in optional_fields:
        builder.set_alignment(
            optional_fields['alignment_method'],
            optional_fields.get('alignment_params')
        )
    if 'filter_type' in optional_fields:
        builder.set_filter(optional_fields['filter_type'])
    if 'modifier' in optional_fields:
        builder.set_modifier(optional_fields['modifier'])
    if 'reference' in optional_fields:
        builder.add_reference(optional_fields['reference'])
    if 'removal' in optional_fields:
        builder.set_removal(optional_fields['removal'])
    if 'reviewer' in optional_fields:
        builder.set_reviewer(optional_fields['reviewer'])
    if 'extras' in optional_fields:
        builder.set_extras(optional_fields['extras'])
    
    # Handle raw optional fields (e.g., alignment='a-maff')
    if 'alignment' in optional_fields and 'alignment_method' not in optional_fields:
        builder._data['alignment'] = optional_fields['alignment']
    if 'filters' in optional_fields and 'filter_type' not in optional_fields:
        builder._data['filters'] = optional_fields['filters']
    if 'modifiers' in optional_fields and 'modifier' not in optional_fields:
        builder._data['modifiers'] = optional_fields['modifiers']
    if 'additions' in optional_fields and 'reference' not in optional_fields:
        builder._data['additions'] = optional_fields['additions']
    if 'removals' in optional_fields and 'removal' not in optional_fields:
        builder._data['removals'] = optional_fields['removals']
    if 'sequencing' in optional_fields:
        builder._data['sequencing'] = optional_fields['sequencing']
    if 'processing' in optional_fields:
        builder._data['processing'] = optional_fields['processing']
    if 'reviewed' in optional_fields and 'reviewer' not in optional_fields:
        builder._data['reviewed'] = optional_fields['reviewed']
    
    filename_obj = builder.build()
    
    # Load sequences
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        try:
            seq_name = parse_sequence_name(record.id)
            sequences.append({
                'record': record,
                'parsed_name': seq_name
            })
        except ValueError as e:
            print(f"Warning: Could not parse sequence name '{record.id}': {e}", file=sys.stderr)
            sequences.append({
                'record': record,
                'parsed_name': None
            })
    
    return {
        'filename': filename_obj,
        'sequences': sequences,
        'original_filename': fasta_path.name
    }


def write_fasta_with_naming(
    filename: FastaFileName,
    sequences: List[Dict[str, any]],
    output_path: Path,
    naming_convention: str = 'standard',
    include_optional: Union[bool, str] = False
) -> Path:
    """
    Write FASTA file with proper naming convention.
    
    Args:
        filename: FastaFileName object
        sequences: List of dicts with 'record' and optionally 'parsed_name'
        output_path: Directory to write file
        naming_convention: Naming convention to use. Supported: 'standard', 'hvtn', 'caprisa'
        include_optional: For 'standard' convention only:
                         - True: Include all available optional fields
                         - False: Include no optional fields (default)
                         - str: Comma-separated flags (e.g., 'a,f,r' for alignment, filters, reviewed)
                               Valid flags: a(alignment), f(filters), m(modifiers), w(additions),
                               o(removals), s(sequencing), p(processing), e(extras), r(reviewed)
    
    Returns:
        Path to written file
    
    Examples:
        # Standard naming convention (obligate fields only)
        write_fasta_with_naming(fn, seqs, Path("."))
        
        # Include all optional fields
        write_fasta_with_naming(fn, seqs, Path("."), include_optional=True)
        
        # Include specific optional fields
        write_fasta_with_naming(fn, seqs, Path("."), include_optional='a,f,r')
        
        # HVTN naming convention
        write_fasta_with_naming(fn, seqs, Path("."), naming_convention='hvtn')
    
    Adding new naming conventions:
        1. Define new convention class inheriting from base naming classes
        2. Implement to_filename() and to_sequence_id() methods
        3. Add convention name to _apply_naming_convention()
        4. Add any convention-specific validation rules
    """
    # Prepare records with proper naming based on convention
    records_to_write = []
    for seq_data in sequences:
        record = seq_data['record']
        parsed_name = seq_data.get('parsed_name')
        
        if parsed_name:
            # Apply naming convention to generate proper ID
            record.id = _apply_naming_convention(parsed_name, naming_convention)
            record.description = ""
        
        records_to_write.append(record)
    
    # Generate output filename based on convention
    output_filename = _generate_filename(filename, naming_convention, include_optional)
    output_file = output_path / output_filename
    
    # Write FASTA file
    SeqIO.write(records_to_write, output_file, "fasta")
    print(f"Wrote {len(records_to_write)} sequences to {output_file} ({naming_convention} convention)")
    
    return output_file


def _apply_naming_convention(seq_name: SequenceName, convention: str) -> str:
    """
    Apply naming convention to generate sequence ID.
    
    Args:
        seq_name: SequenceName object
        convention: Naming convention ('standard', 'hvtn', 'caprisa')
    
    Returns:
        Formatted sequence ID string
    """
    if convention == 'standard':
        # Default convention: Network-Cohort_PTID_Visit_Region_Method-SeqID
        return seq_name.to_sequence_id()
    
    elif convention == 'hvtn':
        # HVTN convention: PTID.Visit.Region.SeqID
        parts = [
            seq_name.ptid,
            seq_name.visit,
            seq_name.region,
        ]
        if seq_name.sequence_id:
            parts.append(seq_name.sequence_id)
        return ".".join(parts)
    
    elif convention == 'caprisa':
        # CAPRISA convention: PTID_Visit_Region_Method_SeqID
        parts = [
            seq_name.ptid,
            seq_name.visit,
            seq_name.region,
            seq_name.sequencing_method
        ]
        if seq_name.sequence_id:
            parts.append(seq_name.sequence_id)
        return "_".join(parts)
    
    else:
        raise ValueError(f"Unknown naming convention: {convention}")


def _generate_filename(filename: FastaFileName, convention: str, include_optional: Union[bool, str] = False) -> str:
    """
    Generate filename based on naming convention.
    
    Args:
        filename: FastaFileName object
        convention: Naming convention ('standard', 'hvtn', 'caprisa')
        include_optional: For 'standard' convention:
                         - True: Include all available optional fields
                         - False: Include no optional fields
                         - str: Comma-separated flags to include
    
    Returns:
        Formatted filename string
    """
    if convention == 'standard':
        # Determine which optional fields to include
        if include_optional is False:
            # Only obligate fields
            parts = [
                f"{filename.network}{filename.protocol}",
                filename.ptid,
                filename.visit,
                filename.region,
                filename.molecule
            ]
        elif include_optional is True:
            # All available optional fields
            return filename.to_filename()
        else:
            # Parse flags string
            flags = [f.strip() for f in include_optional.split(',')]
            
            # Start with obligate fields
            parts = [
                f"{filename.network}{filename.protocol}",
                filename.ptid,
                filename.visit,
                filename.region,
                filename.molecule
            ]
            
            # Add requested optional fields in standard order
            for flag in ['a', 'f', 'm', 'w', 'o', 's', 'p', 'e', 'r']:
                if flag in flags:
                    field_name = OPTIONAL_FIELD_FLAGS[flag]
                    field_value = getattr(filename, field_name, None)
                    if field_value:
                        parts.append(_format_optional_field(flag, field_value))
        
        return "_".join(parts) + ".fasta"
    
    elif convention == 'hvtn':
        # HVTN convention: Protocol_PTID_Visit_Region_Molecule.fasta
        return f"{filename.protocol}_{filename.ptid}_{filename.visit}_{filename.region}_{filename.molecule}.fasta"
    
    elif convention == 'caprisa':
        # CAPRISA convention: PTID_Visit_Region_Molecule.fasta
        return f"{filename.ptid}_{filename.visit}_{filename.region}_{filename.molecule}.fasta"
    
    else:
        raise ValueError(f"Unknown naming convention: {convention}")


# ============================================================================
# Validation and Utility Functions
# ============================================================================

def validate_fasta_file(fasta_path: Path, verbose: bool = True) -> List[str]:
    """
    Validate FASTA file against naming convention.
    
    Args:
        fasta_path: Path to FASTA file
        verbose: Print validation results
        
    Returns:
        List of validation error messages (empty if valid)
    """
    errors = []
    
    # Validate filename
    try:
        filename_obj = parse_fasta_filename(fasta_path.name)
        if verbose:
            print(f"✓ Filename valid: {fasta_path.name}")
    except ValueError as e:
        errors.append(f"Filename error: {e}")
        if verbose:
            print(f"✗ {errors[-1]}")
        return errors
    
    # Validate sequences
    seq_count = 0
    for i, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
        seq_count += 1
        try:
            seq_name = parse_sequence_name(record.id)
            
            # Check consistency between filename and sequence names
            if seq_name.network != filename_obj.network:
                errors.append(f"Sequence {i+1}: network mismatch (file: {filename_obj.network}, seq: {seq_name.network})")
            
            if seq_name.ptid != filename_obj.ptid:
                errors.append(f"Sequence {i+1}: PTID mismatch (file: {filename_obj.ptid}, seq: {seq_name.ptid})")
                
        except ValueError as e:
            errors.append(f"Sequence {i+1} ('{record.id}'): {e}")
    
    if verbose:
        if errors:
            print(f"✗ Found {len(errors)} error(s) in {seq_count} sequences")
            for err in errors:
                print(f"  - {err}")
        else:
            print(f"✓ All {seq_count} sequences valid")
    
    return errors


def convert_naming_convention(
    input_fasta: Path,
    output_dir: Path,
    new_filename: Optional[FastaFileName] = None,
    sequence_transformer: Optional[callable] = None,
    naming_convention: str = 'standard',
    include_optional: Union[bool, str] = False
) -> Path:
    """
    Convert FASTA file to new naming convention.
    
    Args:
        input_fasta: Input FASTA file
        output_dir: Output directory
        new_filename: New filename object (None to keep existing)
        sequence_transformer: Function to transform SequenceName objects
        naming_convention: Target naming convention ('standard', 'hvtn', 'caprisa')
        include_optional: For 'standard' convention: True (all), False (none), or 
                         comma-separated flags (e.g., 'a,f,r')
        
    Returns:
        Path to output file
    """
    data = read_fasta_with_naming(input_fasta)
    
    filename = new_filename if new_filename else data['filename']
    
    if sequence_transformer:
        for seq_data in data['sequences']:
            if seq_data['parsed_name']:
                seq_data['parsed_name'] = sequence_transformer(seq_data['parsed_name'])
    
    output_file = write_fasta_with_naming(
        filename, data['sequences'], output_dir, 
        naming_convention=naming_convention,
        include_optional=include_optional
    )
    
    return output_file


# ============================================================================
# Command Line Interface
# ============================================================================

def main():
    """Command line interface."""
    if len(sys.argv) < 2:
        print("Usage: python fasta_naming.py <command> [arguments]")
        print("\nCommands:")
        print("  validate <file.fasta>              - Validate FASTA file naming")
        print("  parse <file.fasta>                 - Parse and display naming components")
        print("  convert <file.fasta> <convention>  - Convert to different naming convention")
        print("  help                               - Show this help message")
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == "help":
        print(__doc__)
        sys.exit(0)
    
    elif command == "validate":
        if len(sys.argv) < 3:
            print("Error: Please provide a FASTA file to validate", file=sys.stderr)
            sys.exit(1)
        
        fasta_path = Path(sys.argv[2])
        if not fasta_path.exists():
            print(f"Error: File not found: {fasta_path}", file=sys.stderr)
            sys.exit(1)
        
        errors = validate_fasta_file(fasta_path, verbose=True)
        sys.exit(0 if not errors else 1)
    
    elif command == "parse":
        if len(sys.argv) < 3:
            print("Error: Please provide a FASTA file to parse", file=sys.stderr)
            sys.exit(1)
        
        fasta_path = Path(sys.argv[2])
        if not fasta_path.exists():
            print(f"Error: File not found: {fasta_path}", file=sys.stderr)
            sys.exit(1)
        
        try:
            data = read_fasta_with_naming(fasta_path)
            
            print("\n=== Filename Components ===")
            fn = data['filename']
            print(f"Network:    {fn.network}")
            print(f"Protocol:   {fn.protocol}")
            print(f"PTID:       {fn.ptid}")
            print(f"Visit:      {fn.visit}")
            print(f"Region:     {fn.region}")
            print(f"Molecule:   {fn.molecule}")
            
            if fn.alignment:
                print(f"Alignment:  {fn.alignment}")
            if fn.filters:
                print(f"Filters:    {fn.filters}")
            if fn.modifiers:
                print(f"Modifiers:  {fn.modifiers}")
            if fn.reviewed:
                print(f"Reviewed:   {fn.reviewed}")
            
            print(f"\n=== Sequences ({len(data['sequences'])}) ===")
            for i, seq_data in enumerate(data['sequences'][:5], 1):  # Show first 5
                sn = seq_data['parsed_name']
                if sn:
                    print(f"\n{i}. {seq_data['record'].id}")
                    print(f"   Network: {sn.network}, Cohort: {sn.cohort}, PTID: {sn.ptid}")
                    print(f"   Visit: {sn.visit}, Region: {sn.region}")
                    print(f"   Sequencing: {sn.sequencing_method}-{sn.sequence_id}")
                else:
                    print(f"\n{i}. {seq_data['record'].id} (parsing failed)")
            
            if len(data['sequences']) > 5:
                print(f"\n... and {len(data['sequences']) - 5} more sequences")
            
        except Exception as e:
            print(f"Error parsing file: {e}", file=sys.stderr)
            sys.exit(1)
    
    elif command == "convert":
        if len(sys.argv) < 4:
            print("Error: Usage: convert <input.fasta> <convention>", file=sys.stderr)
            print("Supported conventions: standard, hvtn, caprisa", file=sys.stderr)
            sys.exit(1)
        
        input_path = Path(sys.argv[2])
        convention = sys.argv[3]
        
        if not input_path.exists():
            print(f"Error: File not found: {input_path}", file=sys.stderr)
            sys.exit(1)
        
        try:
            output_file = convert_naming_convention(
                input_path,
                input_path.parent,
                naming_convention=convention
            )
            print(f"Converted to {convention} convention: {output_file}")
        except Exception as e:
            print(f"Error converting file: {e}", file=sys.stderr)
            sys.exit(1)
    
    else:
        print(f"Error: Unknown command '{command}'", file=sys.stderr)
        print("Run 'python fasta_naming.py help' for usage", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
