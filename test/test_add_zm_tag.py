import pytest
from src.p_hacking.add_zm_tag import derive_pu_from_qname, extract_zm_from_qname

def test_derive_pu_from_qname():
    # Test cases: (input_query_name, expected_pu)
    test_cases = [
        ("m54006_160901_123456/123456/0_1000", "m54006_160901_123456"),
        ("m12345_678901_234567/1/0_500", "m12345_678901_234567"),
        ("read_without_slashes", "unknown"),
        ("m98765_432109_876543/2/", "m98765_432109_876543"),
        ("/0_100", "unknown"),
        ("", "unknown"),
    ]

    for query_name, expected_pu in test_cases:
        assert derive_pu_from_qname(query_name) == expected_pu
        
def test_extract_zm_from_qname():
    # Test cases: (input_query_name, expected_zm)
    test_cases = [
        ("movie/42", 42),
        ("sample/0", 0),
        ("test/9999", 9999),
        ("no_slash", None),
        ("invalid/abc", None),
        ("another_invalid/", None),
        ("/123", 123),
        ("", None),
    ]

    for query_name, expected_zm in test_cases:
        assert extract_zm_from_qname(query_name) == expected_zm