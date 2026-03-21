"""Shared path helpers for the analysis pipeline.

All paths are anchored at the repository root, regardless of where
scripts are executed from.
"""
import os

# src/ directory where this file lives
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# Repository root (one level up from src/)
ROOT_DIR = os.path.abspath(os.path.join(BASE_DIR, ".."))

RAW_DIR     = os.path.join(ROOT_DIR, "data", "raw")
DERIVED_DIR = os.path.join(ROOT_DIR, "data", "derived")
RESULTS_DIR = os.path.join(ROOT_DIR, "results")


def raw_path(*parts):
    return os.path.join(RAW_DIR, *parts)

def derived_path(*parts):
    return os.path.join(DERIVED_DIR, *parts)

def results_path(*parts):
    return os.path.join(RESULTS_DIR, *parts)
