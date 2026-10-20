import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import fetch_interactions as fi  # noqa: E402


class FakeWriter:
    def __init__(self):
        self.rows = []

    def write_row(self, **row):
        self.rows.append(row)


def row(uid, gene, partner):
    return {"uniqueId": uid, "geneId": gene, "partnerGeneId": partner}


def run(batches):
    writer, seen, pending = FakeWriter(), set(), {}
    for batch in batches:
        fi._write_new(writer, batch, seen, pending)
    fi._flush_pending(writer, pending, seen)
    return [(r["geneId"], r["partnerGeneId"]) for r in writer.rows]


def test_direction_does_not_depend_on_arrival_order():
    a_first = run([[row("u1", "A", "B")], [row("u1", "B", "A")]])
    b_first = run([[row("u1", "B", "A")], [row("u1", "A", "B")]])
    assert a_first == b_first == [("A", "B")]


def test_single_non_canonical_copy_is_kept_after_flush():
    assert run([[row("u1", "B", "A")]]) == [("B", "A")]


def test_flush_is_sorted_by_unique_id():
    out = run([[row("u2", "D", "C"), row("u1", "B", "A")]])
    assert out == [("B", "A"), ("D", "C")]


def test_rows_without_unique_id_are_written():
    assert run([[row("", "B", "A"), row("", "B", "A")]]) == [("B", "A"), ("B", "A")]
