# DRAFT — for ecere/dggal (delete this file after filing)

Title:

> ISEA7H: `getZoneFromWGS84Centroid` returns the null zone for some points at level 15

---

When using ISEA7H, `getZoneFromWGS84Centroid` returns the null zone
(`0xffffffffffffffff`) for some points at level 15, even though the same
points resolve to valid zones at every other level.

Here's one example:

```python
import dggal as al
al.pydggal_setup(al.Application())
isea = al.ISEA7H()

p = al.GeoPoint(63.466273, 11.200161)
z = isea.getZoneFromWGS84Centroid(15, p)

print(hex(z))
print(isea.getZoneTextID(z))
for v in isea.getZoneWGS84Vertices(z):
    print((float(v.lat), float(v.lon)))
```

Output:

```
0xffffffffffffffff
(null)
(0.0, 0.0)
(0.0, 0.0)
(0.0, 0.0)
(0.0, 0.0)
(0.0, 0.0)
(0.0, 0.0)
```

The same point resolves to a valid zone at levels 0–14 and 16–19, and at
every level of IVEA7H and rHEALPix. A nearby point works fine at level 15
(`(63.476273, 11.200161)` → `P8-D2B25A692-F`), so this seems specific to
level 15 rather than a general face-boundary issue — though the failing
point is near an icosahedron edge.

Found via uniform random sampling; very roughly one point per million draws
at level 15 seems to hit this. A secondary papercut: zone APIs called on the
null zone return plausible-looking values (`(null)`, `(0, 0)` vertices)
rather than erroring, which makes the failure easy to ingest silently.

Environment: dggal 0.0.6, Python 3.13.12; reproduced on Linux (x86_64) and
macOS 26.5.1 (x86_64 wheel via Rosetta on arm64).
