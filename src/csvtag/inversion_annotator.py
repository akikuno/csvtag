import cstag
import midsv


midsv.transform(sam, midsv=False, cssplit=True, qscore=False, keep=set(["FLAG"]))
cstag.read_sam

def annotate_inversion(cs_tags: str | list[str], positions: int | list[int], flags: int | list[int]) -> str:
    records = []
    for cs_tag, pos, flag in zip(cs_tags, positions, flags):
        records.append({"cs_tag": cs_tag, "pos": pos, "flag": flag})

    records.sort(key=lambda x: x["pos"])
