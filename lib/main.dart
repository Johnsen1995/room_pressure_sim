// main.dart — Room Pressure Simulator (Flutter)
//
// Interactions
// - ROOM: drag corners. Double‑tap an edge to add a corner. Long‑press a corner to delete.
// - SPEAKERS: tap to add. Long‑press near a speaker to remove.
// - SEAT: tap to place the seat.
// - AREA: drag to draw allowed speaker area. Long‑press inside the area to clear.
//
// Simulation model (very simplified):
// - Uses "image sources" (mirror reflections) of each speaker across each wall segment.
// - For reflection order N, it generates sources by repeatedly reflecting across all wall edges.
// - Wall reflection factor controls amplitude per order (wallRefl^order).
// - SPL is computed as magnitude of summed complex waves at each grid point.
//
// Heatmap colors are RELATIVE to the loudest point in the current simulation:
// green ≈ -30 dB, yellow ≈ -15 dB, red ≈ 0 dB (max).

import 'dart:async';
import 'dart:convert';
import 'dart:io';
import 'dart:math' as math;
import 'dart:typed_data';
import 'dart:ui' as ui;

import 'package:file_picker/file_picker.dart';
import 'package:flutter/material.dart';
import 'package:flutter/services.dart';
import 'package:path_provider/path_provider.dart';

const double kHeatmapSpanDb =
    30.0; // Must match simulation color normalization (-30..0 dB relative to max)

void main() {
  WidgetsFlutterBinding.ensureInitialized();
  runApp(const RoomSimBootstrap());
}

class RoomSimApp extends StatelessWidget {
  const RoomSimApp({super.key});

  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      title: 'Room Pressure Simulator',
      debugShowCheckedModeBanner: false,
      theme: ThemeData(useMaterial3: true),
      home: const SimulatorPage(),
    );
  }
}

class RoomSimBootstrap extends StatefulWidget {
  const RoomSimBootstrap({super.key});

  @override
  State<RoomSimBootstrap> createState() => _RoomSimBootstrapState();
}

class _RoomSimBootstrapState extends State<RoomSimBootstrap> {
  late Future<bool> _seenFuture;

  @override
  void initState() {
    super.initState();
    _seenFuture = _IntroSeenStore.isSeen();
  }

  @override
  Widget build(BuildContext context) {
    return FutureBuilder<bool>(
      future: _seenFuture,
      builder: (context, snap) {
        final bool seen = snap.data == true;
        return MaterialApp(
          title: 'Room Pressure Simulator',
          debugShowCheckedModeBanner: false,
          theme: ThemeData(useMaterial3: true),
          home: seen
              ? const SimulatorPage()
              : IntroScreen(
                  firstLaunch: true,
                  onDone: () async {
                    await _IntroSeenStore.markSeen();
                    if (!mounted) return;
                    setState(() {
                      _seenFuture = Future<bool>.value(true);
                    });
                  },
                ),
        );
      },
    );
  }
}

/// Tiny persistent flag to show the intro once on first launch.
/// Uses a JSON file in the application documents directory to avoid extra deps.
class _IntroSeenStore {
  static const String _fileName = 'intro_seen.json';

  static Future<File> _file() async {
    final dir = await getApplicationDocumentsDirectory();
    return File('${dir.path}/$_fileName');
  }

  static Future<bool> isSeen() async {
    try {
      final f = await _file();
      if (!await f.exists()) return false;
      final raw = await f.readAsString();
      final obj = jsonDecode(raw);
      return obj is Map && obj['seen'] == true;
    } catch (_) {
      // Fail open: don't block the app if storage is unavailable.
      return false;
    }
  }

  static Future<void> markSeen() async {
    try {
      final f = await _file();
      await f.writeAsString(jsonEncode({'seen': true}));
    } catch (_) {
      // Ignore.
    }
  }
}

class IntroScreen extends StatelessWidget {
  const IntroScreen({
    super.key,
    required this.firstLaunch,
    required this.onDone,
  });

  final bool firstLaunch;
  final Future<void> Function() onDone;

  @override
  Widget build(BuildContext context) {
    final String primaryLabel = firstLaunch ? 'Got it' : 'Close';

    return Scaffold(
      appBar: AppBar(
        title: const Text('Room Pressure Simulator'),
        actions: [
          if (!firstLaunch)
            IconButton(
              tooltip: 'Close',
              icon: const Icon(Icons.close),
              onPressed: () => Navigator.of(context).maybePop(),
            ),
        ],
      ),
      body: SafeArea(
        child: Padding(
          padding: const EdgeInsets.fromLTRB(20, 16, 20, 20),
          child: Column(
            crossAxisAlignment: CrossAxisAlignment.start,
            children: [
              const Text(
                'What this shows',
                style: TextStyle(fontSize: 16, fontWeight: FontWeight.w800),
              ),
              const SizedBox(height: 8),
              const Text(
                'This app visualizes how speaker placement affects relative sound pressure in a room.\n\n'
                'The heatmap is normalized to a reference level (0 dB), chosen in Settings.\n'
                'Use it to compare placements and build intuition — not to predict exact measurements.',
                style: TextStyle(fontSize: 13, height: 1.35),
              ),
              const SizedBox(height: 16),
              const Text(
                'What it does not model',
                style: TextStyle(fontSize: 16, fontWeight: FontWeight.w800),
              ),
              const SizedBox(height: 8),
              const Text(
                '• room reflections or standing waves\n'
                '• wall absorption or room treatment\n'
                '• detailed speaker or listener directivity\n\n'
                'Complex room shapes are shown visually; acoustic blocking is approximate.',
                style: TextStyle(fontSize: 13, height: 1.35),
              ),
              const SizedBox(height: 16),
              const Text(
                'About suggestions',
                style: TextStyle(fontSize: 16, fontWeight: FontWeight.w800),
              ),
              const SizedBox(height: 8),
              const Text(
                'Suggested/optimized placement (when enabled) is biased toward realistic setups '
                '(e.g. near walls, away from seating). Manual placement is always unrestricted.',
                style: TextStyle(fontSize: 13, height: 1.35),
              ),
              const Spacer(),
              SizedBox(
                width: double.infinity,
                child: FilledButton(
                  onPressed: () async {
                    if (firstLaunch) {
                      await onDone();
                    } else {
                      Navigator.of(context).maybePop();
                    }
                  },
                  child: Text(primaryLabel),
                ),
              ),
            ],
          ),
        ),
      ),
    );
  }
}

enum StereoPreset { typical, wide, tooWide }

enum PlacementStyle { realistic, free }

enum EditMode { room, speakers, seat, area }

enum SeatType { couch, chair }

class SimulatorPage extends StatefulWidget {
  const SimulatorPage({super.key});

  @override
  State<SimulatorPage> createState() => _SimulatorPageState();
}

class _SimulatorPageState extends State<SimulatorPage> {
  // ---------- Model ----------
  List<Offset> poly = const [
    Offset(0, 0),
    Offset(5, 0),
    Offset(5, 4),
    Offset(0, 4),
  ].toList();

  List<Offset> speakers = [];
  Offset? seat;
  SeatType seatType = SeatType.couch;

  Rect? allowedArea;

  // Area drag state
  Offset? _areaDragStartWorld;

  // ---------- Modes ----------
  EditMode mode = EditMode.room;
  int? draggingVertex;

  // Prevent tap-up action after a long-press action
  bool _ignoreNextTapUp = false;

  // ---------- Simulation params ----------
  double freqHz = 65.0;
  String optFreqs = "35,45,55,65";
  double gain = 1.0;
  double wallRefl = 0.85;
  int reflections = 2;
  double gridStep = 0.12;

  // ---------- Optimization params ----------
  int optSpeakers = 2;
  int optIters = 800;

  // ---------- Suggestion style (Option B) ----------
  PlacementStyle placementStyle = PlacementStyle.realistic;

  // Heatmap normalization: strict max vs robust percentile (prevents tiny hotspots dominating).
  bool robustHeatmapScale = false;

  // Cached sample points for optimizer (to estimate room max without full grid)
  List<Offset>? _optSamples;

  // ---------- Presets / Defaults ----------
  void _applyPreset(StereoPreset preset) {
    if (poly.length < 4) return;

    // Assume the first edge (0->1) is the "front" wall for presets.
    final p0 = poly[0];
    final p1 = poly[1];
    final p3 = poly.length > 3 ? poly[3] : poly.last;

    final xAxisRaw = p1 - p0;
    final yAxisRaw = p3 - p0;
    final width = xAxisRaw.distance;
    final depth = yAxisRaw.distance;

    if (width < 1e-6 || depth < 1e-6) return;

    final xAxis = xAxisRaw / width;
    final yAxis = yAxisRaw / depth;

    const inset = 0.35; // meters from walls (keeps placements "sane")
    final fxInset = (inset / width).clamp(0.02, 0.15);
    final fyInset = (inset / depth).clamp(0.02, 0.20);

    double leftFx;
    double rightFx;
    switch (preset) {
      case StereoPreset.typical:
        leftFx = 0.25;
        rightFx = 0.75;
        break;
      case StereoPreset.wide:
        leftFx = 0.18;
        rightFx = 0.82;
        break;
      case StereoPreset.tooWide:
        leftFx = 0.05;
        rightFx = 0.95;
        break;
    }

    // Clamp away from corners a bit (even for "too wide", keep it inside).
    leftFx = leftFx.clamp(fxInset, 1.0 - fxInset);
    rightFx = rightFx.clamp(fxInset, 1.0 - fxInset);

    Offset pt(double fx, double fy) =>
        p0 + xAxis * (fx * width) + yAxis * (fy * depth);

    final spkY = fyInset;
    final seatY = (0.75).clamp(fyInset * 2, 1.0 - fyInset);

    setState(() {
      speakers = [pt(leftFx, spkY), pt(rightFx, spkY)];
      seat = pt(0.5, seatY);

      // Presets are a "layout change" -> clear derived overlays/readouts.
      _heatImage?.dispose();
      _heatImage = null;
      _lastMaxSpl = null;
      _lastSeatRelDb = null;
    });
  }

  void _resetToDefaults() {
    setState(() {
      poly = const [
        Offset(0, 0),
        Offset(5, 0),
        Offset(5, 4),
        Offset(0, 4),
      ].toList();

      allowedArea = null;
      mode = EditMode.room;
      seatType = SeatType.couch;

      _heatImage?.dispose();
      _heatImage = null;
      _lastMaxSpl = null;
      _lastSeatRelDb = null;
    });

    // Apply a "sane" starting layout.
    _applyPreset(StereoPreset.typical);
  }

  // ---------- Rendering ----------
  Rect _viewWorld = const Rect.fromLTWH(-0.5, -0.5, 6.0, 5.0);
  Rect _heatWorld = const Rect.fromLTWH(0, 0, 5, 4);
  ui.Image? _heatImage;
  bool _busy = false;
  double? _lastMaxSpl;

  double? _lastSeatRelDb; // Relative dB (0 dB at selected heatmap reference)

  // UI
  bool _showLegendDetails = true;
  bool _occlusionEnabled = false;
  bool _softDiffraction = true;
  // ---------- Tutorial ----------
  bool showTutorial = false;
  int tutorialStep = 0;

  final GlobalKey _kModeRow = GlobalKey();
  final GlobalKey _kCanvas = GlobalKey();
  final GlobalKey _kSimBtn = GlobalKey();
  final GlobalKey _kOptBtn = GlobalKey();
  final GlobalKey _kSettingsBtn = GlobalKey();

  static const double _pickRadiusPx = 18.0;
  static const double _minEdgeInsertDistPx = 22.0;

  @override
  void initState() {
    super.initState();
    _viewWorld = _polyBounds(poly).inflate(0.7);
    _heatWorld = _polyBounds(poly).inflate(0.001);
    _occlusionEnabled = !_isConvex(poly);
  }

  @override
  void dispose() {
    _heatImage?.dispose();
    super.dispose();
  }

  // ---------- Geometry ----------
  Rect _polyBounds(List<Offset> p) {
    double minX = p.first.dx, maxX = p.first.dx;
    double minY = p.first.dy, maxY = p.first.dy;
    for (final v in p) {
      minX = math.min(minX, v.dx);
      maxX = math.max(maxX, v.dx);
      minY = math.min(minY, v.dy);
      maxY = math.max(maxY, v.dy);
    }
    return Rect.fromLTRB(minX, minY, maxX, maxY);
  }

  bool _pointInPoly(Offset p, List<Offset> poly) {
    bool inside = false;
    for (int i = 0; i < poly.length; i++) {
      final a = poly[i];
      final b = poly[(i + 1) % poly.length];
      final cond =
          ((a.dy > p.dy) != (b.dy > p.dy)) &&
          (p.dx <
              (b.dx - a.dx) * (p.dy - a.dy) / ((b.dy - a.dy) + 1e-15) + a.dx);
      if (cond) inside = !inside;
    }
    return inside;
  }

  double _cross(Offset a, Offset b) => a.dx * b.dy - a.dy * b.dx;

  bool _isConvex(List<Offset> p) {
    if (p.length < 4) return true; // triangles always convex
    double? sign;
    for (int i = 0; i < p.length; i++) {
      final a = p[i];
      final b = p[(i + 1) % p.length];
      final c = p[(i + 2) % p.length];
      final ab = b - a;
      final bc = c - b;
      final z = _cross(ab, bc);
      if (z.abs() < 1e-10) continue;
      final s = z.sign;
      sign ??= s;
      if (sign != s) return false;
    }
    return true;
  }

  double _distPointToSegment(Offset p, Offset a, Offset b) {
    final ab = b - a;
    final ab2 = ab.distanceSquared;
    if (ab2 < 1e-12) return (p - a).distance;
    final ap = p - a;
    double t = (ap.dx * ab.dx + ap.dy * ab.dy) / ab2;
    t = t.clamp(0.0, 1.0);
    final proj = a + ab * t;
    return (p - proj).distance;
  }

  bool _segmentsProperlyIntersect(Offset p1, Offset p2, Offset q1, Offset q2) {
    // Proper intersection (not merely touching at endpoints).
    // Uses orientation tests with epsilon.
    const eps = 1e-10;

    double orient(Offset a, Offset b, Offset c) {
      return (b.dx - a.dx) * (c.dy - a.dy) - (b.dy - a.dy) * (c.dx - a.dx);
    }

    bool onSegment(Offset a, Offset b, Offset c) {
      // c on segment ab (colinear case)
      return math.min(a.dx, b.dx) - eps <= c.dx &&
          c.dx <= math.max(a.dx, b.dx) + eps &&
          math.min(a.dy, b.dy) - eps <= c.dy &&
          c.dy <= math.max(a.dy, b.dy) + eps;
    }

    final o1 = orient(p1, p2, q1);
    final o2 = orient(p1, p2, q2);
    final o3 = orient(q1, q2, p1);
    final o4 = orient(q1, q2, p2);

    // General case
    if ((o1 > eps && o2 < -eps) || (o1 < -eps && o2 > eps)) {
      if ((o3 > eps && o4 < -eps) || (o3 < -eps && o4 > eps)) return true;
    }

    // Colinear / touching cases: treat as not "proper" if it's just endpoint touch
    if (o1.abs() <= eps && onSegment(p1, p2, q1)) return true;
    if (o2.abs() <= eps && onSegment(p1, p2, q2)) return true;
    if (o3.abs() <= eps && onSegment(q1, q2, p1)) return true;
    if (o4.abs() <= eps && onSegment(q1, q2, p2)) return true;

    return false;
  }

  bool _segmentStaysInsidePoly(Offset a, Offset b, List<Offset> p) {
    // Assumes endpoints are inside polygon.
    // If the segment crosses any polygon edge at a non-endpoint location, it leaves the polygon.
    const eps = 1e-8;

    for (int i = 0; i < p.length; i++) {
      final e1 = p[i];
      final e2 = p[(i + 1) % p.length];

      // Skip edges that share an endpoint with the segment (allow touching at vertices).
      if ((e1 - a).distance < eps ||
          (e2 - a).distance < eps ||
          (e1 - b).distance < eps ||
          (e2 - b).distance < eps) {
        continue;
      }

      if (_segmentsProperlyIntersect(a, b, e1, e2)) return false;
    }
    return true;
  }

  double _occlusionPenaltyDb(Offset a, Offset b) {
    // Heuristic "diffraction": if the blocked path passes near a corner, reduce penalty.
    // Near corner -> -8 dB, far from corners -> -18 dB.
    const near = 0.20; // m
    const far = 1.00; // m
    const nearDb = -8.0;
    const farDb = -18.0;

    double minD = double.infinity;
    for (final v in poly) {
      minD = math.min(minD, _distPointToSegment(v, a, b));
    }

    final t = ((minD - near) / (far - near)).clamp(0.0, 1.0);
    return nearDb + (farDb - nearDb) * t;
  }

  Offset _reflectAcrossLine(Offset p, Offset a, Offset b) {
    final vx = b.dx - a.dx;
    final vy = b.dy - a.dy;
    final vv = vx * vx + vy * vy;
    if (vv.abs() < 1e-12) return p;
    final t = ((p.dx - a.dx) * vx + (p.dy - a.dy) * vy) / vv;
    final proj = Offset(a.dx + t * vx, a.dy + t * vy);
    return Offset(2 * proj.dx - p.dx, 2 * proj.dy - p.dy);
  }

  double _log10(double x) => math.log(x) / math.ln10;

  // ---------- World/screen transforms ----------
  Offset _worldToScreen(Offset w, Size size, Rect view) {
    final sx = (w.dx - view.left) / view.width * size.width;
    final sy = (w.dy - view.top) / view.height * size.height;
    return Offset(sx, sy);
  }

  Offset _screenToWorld(Offset s, Size size, Rect view) {
    final x = view.left + (s.dx / size.width) * view.width;
    final y = view.top + (s.dy / size.height) * view.height;
    return Offset(x, y);
  }

  void _refreshView() {
    _viewWorld = _polyBounds(poly).inflate(0.7);
  }

  void _refreshViewSetState() {
    setState(() => _viewWorld = _polyBounds(poly).inflate(0.7));
  }

  // ---------- Room editing helpers ----------
  int? _pickVertex(Offset screenPos, Size size) {
    int? best;
    double bestD = double.infinity;
    for (int i = 0; i < poly.length; i++) {
      final s = _worldToScreen(poly[i], size, _viewWorld);
      final d2 = (s - screenPos).distanceSquared;
      if (d2 < bestD) {
        bestD = d2;
        best = i;
      }
    }
    if (best == null) return null;
    return (math.sqrt(bestD) <= _pickRadiusPx) ? best : null;
  }

  int? _pickEdgeForInsert(Offset screenPos, Size size) {
    int? bestEdge;
    double bestDist = double.infinity;

    for (int i = 0; i < poly.length; i++) {
      final aW = poly[i];
      final bW = poly[(i + 1) % poly.length];
      final a = _worldToScreen(aW, size, _viewWorld);
      final b = _worldToScreen(bW, size, _viewWorld);

      final ap = screenPos - a;
      final ab = b - a;
      final ab2 = ab.distanceSquared;
      if (ab2 < 1e-9) continue;
      double t = (ap.dx * ab.dx + ap.dy * ab.dy) / ab2;
      t = t.clamp(0.0, 1.0);
      final proj = Offset(a.dx + t * ab.dx, a.dy + t * ab.dy);
      final dist = (proj - screenPos).distance;

      if (dist < bestDist) {
        bestDist = dist;
        bestEdge = i;
      }
    }

    if (bestEdge == null) return null;
    return (bestDist <= _minEdgeInsertDistPx) ? bestEdge : null;
  }

  void _insertCornerAtEdge(int edgeIndex, Offset worldPoint) {
    setState(() {
      poly = [
        ...poly.sublist(0, edgeIndex + 1),
        worldPoint,
        ...poly.sublist(edgeIndex + 1),
      ];
    });
  }

  void _deleteCorner(int idx) {
    if (poly.length <= 3) return;
    setState(() => poly.removeAt(idx));
  }

  // ---------- Image sources / SPL ----------
  List<_ImageSource> _buildImageSources(Offset speaker) {
    final edges = <(Offset, Offset)>[];
    for (int i = 0; i < poly.length; i++) {
      edges.add((poly[i], poly[(i + 1) % poly.length]));
    }

    final sources = <_ImageSource>[_ImageSource(speaker, 0, 1.0)];
    var current = <_ImageSource>[_ImageSource(speaker, 0, 1.0)];

    for (int order = 1; order <= reflections; order++) {
      final next = <_ImageSource>[];
      for (final src in current) {
        for (final e in edges) {
          final r = _reflectAcrossLine(src.pos, e.$1, e.$2);
          next.add(
            _ImageSource(r, order, math.pow(wallRefl, order).toDouble()),
          );
        }
      }
      sources.addAll(next);
      current = next;
    }
    return sources;
  }

  double _splAtPoint({
    required List<Offset> spks,
    required Offset p,
    required double freq,
  }) {
    const c = 343.0;
    final k = 2 * math.pi * freq / c;
    const eps = 1e-6;

    double sumRe = 0.0;
    double sumIm = 0.0;

    for (final spk in spks) {
      final bool applyOcc = _occlusionEnabled;
      final bool losClear = !applyOcc || _segmentStaysInsidePoly(spk, p, poly);
      final double occDb = (!losClear)
          ? (_softDiffraction ? _occlusionPenaltyDb(spk, p) : -18.0)
          : 0.0;
      final double occAmp = (!losClear)
          ? math.pow(10.0, occDb / 20.0).toDouble()
          : 1.0;

      final imgs = _buildImageSources(spk);
      for (final src in imgs) {
        final dx = p.dx - src.pos.dx;
        final dy = p.dy - src.pos.dy;
        final r = math.sqrt(dx * dx + dy * dy) + eps;
        final baseAmp = gain * src.amp / math.sqrt(r);
        final amp = (src.order == 0) ? (baseAmp * occAmp) : baseAmp;
        final phase = -k * r;
        sumRe += amp * math.cos(phase);
        sumIm += amp * math.sin(phase);
      }
    }

    final mag = math.sqrt(sumRe * sumRe + sumIm * sumIm);
    return 20.0 * _log10(mag + 1e-12);
  }

  List<double> _parseFreqs(String s, double fallback) {
    final parts = s.split(',');
    final out = <double>[];
    for (final p in parts) {
      final v = double.tryParse(p.trim());
      if (v != null && v > 0) out.add(v);
    }
    return out.isEmpty ? [fallback] : out;
  }

  // ---------- Simulation ----------
  Future<void> _simulate() async {
    if (speakers.isEmpty) {
      _toast("Add at least 1 speaker.");
      return;
    }

    setState(() => _busy = true);
    await Future<void>.delayed(const Duration(milliseconds: 10));

    try {
      final b = _polyBounds(poly).inflate(0.001);
      _heatWorld = b;

      final w = ((b.width / gridStep).ceil() + 1).clamp(30, 450);
      final h = ((b.height / gridStep).ceil() + 1).clamp(30, 450);

      final z = Float32List(w * h);
      double zMax = -1e18;

      for (int j = 0; j < h; j++) {
        final y = b.top + j * (b.height / (h - 1));
        for (int i = 0; i < w; i++) {
          final x = b.left + i * (b.width / (w - 1));
          final p = Offset(x, y);
          final idx = j * w + i;

          if (!_pointInPoly(p, poly)) {
            z[idx] = floatNaN;
            continue;
          }
          final spl = _splAtPoint(spks: speakers, p: p, freq: freqHz);
          z[idx] = spl;
          if (spl > zMax) zMax = spl;
        }
      }

      // Choose normalization reference for the heatmap.
      // Strict: absolute max. Robust: high percentile (prevents tiny hotspots dominating the scale).
      double zRef = zMax;
      if (robustHeatmapScale) {
        final values = <double>[];
        // Down-sample for speed on dense grids.
        final int stride = (z.length / 5000).ceil().clamp(1, 999999);
        for (int k = 0; k < z.length; k += stride) {
          final v = z[k];
          if (!v.isNaN && v.isFinite) values.add(v);
        }
        if (values.length >= 20) {
          values.sort();
          final idx = ((values.length - 1) * 0.98).round().clamp(
            0,
            values.length - 1,
          );
          zRef = values[idx];
          if (!zRef.isFinite) zRef = zMax;
        }
      }

      // Seat readout (relative dB) should match heatmap normalization.
      double? seatRelDb;
      if (seat != null && _pointInPoly(seat!, poly)) {
        final seatSpl = _splAtPoint(spks: speakers, p: seat!, freq: freqHz);
        seatRelDb = seatSpl - zRef; // 0 dB at the selected heatmap reference.
      }

      final rgba = Uint8List(w * h * 4);
      for (int idx = 0; idx < w * h; idx++) {
        final o = idx * 4;
        final spl = z[idx];
        if (spl.isNaN) {
          rgba[o + 0] = 0;
          rgba[o + 1] = 0;
          rgba[o + 2] = 0;
          rgba[o + 3] = 0;
          continue;
        }
        final rel = (spl - zRef).clamp(-kHeatmapSpanDb, 0.0);
        final t = (rel + kHeatmapSpanDb) / kHeatmapSpanDb;

        int r, g, bC;
        if (t < 0.5) {
          final u = t / 0.5;
          r = (255 * u).round();
          g = (176 + (79 * u)).round();
          bC = 0;
        } else {
          final u = (t - 0.5) / 0.5;
          r = 255;
          g = (255 * (1 - u)).round();
          bC = 0;
        }

        rgba[o + 0] = r;
        rgba[o + 1] = g;
        rgba[o + 2] = bC;
        rgba[o + 3] = 255;
      }

      final completer = Completer<ui.Image>();
      ui.decodeImageFromPixels(
        rgba,
        w,
        h,
        ui.PixelFormat.rgba8888,
        (img) => completer.complete(img),
      );
      final img = await completer.future;

      _heatImage?.dispose();
      setState(() {
        _heatImage = img;
        _lastMaxSpl = zRef;
        _lastSeatRelDb = seatRelDb;
      });
    } catch (e, st) {
      debugPrint('Simulate failed: $e');
      debugPrintStack(stackTrace: st);
      _toast('Simulation failed (see terminal).');
    } finally {
      if (mounted) setState(() => _busy = false);
    }
  }

  // ---------- Optimization ----------
  Offset? _sampleInPoly({required int tries, Rect? within}) {
    final rng = math.Random();
    final b = within ?? _polyBounds(poly);
    for (int k = 0; k < tries; k++) {
      final x = b.left + rng.nextDouble() * b.width;
      final y = b.top + rng.nextDouble() * b.height;
      final p = Offset(x, y);
      if (_pointInPoly(p, poly)) return p;
    }
    return null;
  }

  double _scoreAtSeat(List<Offset> cand) {
    // Optimizer score: prioritize maximizing the seat level relative to the room max
    // (i.e. make the seat as close to 0 dB on the heatmap as possible).
    if (seat == null) return double.infinity;

    final freqs = _parseFreqs(optFreqs, freqHz);

    // Build a small, stable set of sample points used to estimate the room maximum.
    _optSamples ??= _buildOptSamples(count: 36);

    final rels = <double>[];
    for (final f in freqs) {
      final seatSpl = _splAtPoint(spks: cand, p: seat!, freq: f);
      final roomMax = _estimateRoomMax(
        spks: cand,
        freq: f,
        samples: _optSamples!,
      );
      rels.add(seatSpl - roomMax); // <= 0.0 ideally close to 0
    }

    final meanRel = rels.reduce((a, b) => a + b) / rels.length;
    double varSum = 0.0;
    for (final v in rels) {
      varSum += (v - meanRel) * (v - meanRel);
    }
    final stdRel = math.sqrt(varSum / rels.length);

    // Seat-first: maximize meanRel (closest to 0). Mild penalty for frequency unevenness.
    return (-meanRel) + 0.20 * stdRel;
  }

  List<Offset> _buildOptSamples({required int count}) {
    final out = <Offset>[];
    for (int i = 0; i < count; i++) {
      final p = _sampleInPoly(tries: 6000, within: null);
      if (p == null) break;
      out.add(p);
    }
    // Fallback: include centroid if sampling fails for any reason.
    if (out.isEmpty) {
      Offset c = Offset.zero;
      for (final p in poly) {
        c += p;
      }
      c = c / poly.length.toDouble();
      out.add(c);
    }
    return out;
  }

  double _estimateRoomMax({
    required List<Offset> spks,
    required double freq,
    required List<Offset> samples,
  }) {
    double maxV = -1e30;

    // Base sampling across the room interior.
    for (final p in samples) {
      final v = _splAtPoint(spks: spks, p: p, freq: freq);
      if (v > maxV) maxV = v;
    }

    // Add a few near-speaker points to better catch hot spots without evaluating a full grid.
    const double r = 0.25; // meters
    final dirs = <Offset>[
      const Offset(1, 0),
      const Offset(-1, 0),
      const Offset(0, 1),
      const Offset(0, -1),
    ];
    for (final s in spks) {
      for (final d in dirs) {
        final q = s + d * r;
        if (_pointInPoly(q, poly)) {
          final v = _splAtPoint(spks: spks, p: q, freq: freq);
          if (v > maxV) maxV = v;
        }
      }
    }

    return maxV;
  }

  // ---------- Option B: realism-biased suggestion scoring ----------
  static const double _seatKeepoutM = 0.55; // optimizer only
  static const double _spkKeepoutM = 0.45; // optimizer only

  double _distToNearestWall(Offset p) {
    double best = double.infinity;
    for (int i = 0; i < poly.length; i++) {
      final a = poly[i];
      final b = poly[(i + 1) % poly.length];
      best = math.min(best, _distPointToSegment(p, a, b));
    }
    return best;
  }

  bool _candidateSane(List<Offset> cand) {
    // Keepout around seat (optimizer only). Realistic mode keeps a wider buffer.
    final double seatKeepout =
        _seatKeepoutM +
        (placementStyle == PlacementStyle.realistic ? 0.45 : 0.0);
    if (seat != null) {
      for (final s in cand) {
        if ((s - seat!).distance < seatKeepout) return false;
      }
    }
    // Keepout between speakers.
    for (int i = 0; i < cand.length; i++) {
      for (int j = i + 1; j < cand.length; j++) {
        if ((cand[i] - cand[j]).distance < _spkKeepoutM) return false;
      }
    }
    return true;
  }

  double _placementPenalty(List<Offset> cand) {
    // Manual placement is always unrestricted. This penalty is ONLY used by Optimize/Suggest.
    if (placementStyle == PlacementStyle.free) return 0.0;

    // Realistic style: strongly prefer wall-adjacent placement (speakers in the middle of the room
    // is usually not realistic), but still keep it a SOFT bias (no hard constraints).
    // Also gently discourage being unrealistically glued to a wall.
    const double dPref = 0.35; // beyond this, penalty increases (meters)
    const double dMin = 0.08; // under this, small penalty (meters)

    double penalty = 0.0;

    for (final p in cand) {
      final d = _distToNearestWall(p);

      if (d > dPref) {
        // Quadratic growth: pushes candidates back toward the walls, without forbidding interior points.
        final x = (d - dPref);
        penalty += 1.20 * x * x;
      }
      if (d < dMin) {
        final x = (dMin - d);
        penalty += 0.20 * x * x;
      }
    }

    return penalty;
  }

  Future<void> _optimize() async {
    if (seat == null) {
      _toast("Place the seat first.");
      return;
    }

    setState(() => _busy = true);
    await Future<void>.delayed(const Duration(milliseconds: 10));

    try {
      _optSamples = null;

      final n = optSpeakers.clamp(1, 6);
      final iters = optIters.clamp(50, 20000);

      List<Offset>? best;
      double bestScore = double.infinity;

      for (int it = 0; it < iters; it++) {
        final cand = <Offset>[];
        bool ok = true;

        for (int k = 0; k < n; k++) {
          final p = _sampleInPoly(tries: 8000, within: allowedArea);
          if (p == null) {
            ok = false;
            break;
          }
          cand.add(p);
        }
        if (!ok) continue;

        if (!_candidateSane(cand)) continue;

        final sc = _scoreAtSeat(cand) + _placementPenalty(cand);
        if (sc < bestScore) {
          bestScore = sc;
          best = cand;
        }
      }

      if (best == null) {
        setState(() => _busy = false);
        _toast(
          "No solution found. Increase iterations or enlarge allowed area.",
        );
        return;
      }

      setState(() {
        speakers = best!;
      });

      await _simulate();
      _toast("Optimized.");
    } catch (e, st) {
      debugPrint('Optimize failed: $e');
      debugPrintStack(stackTrace: st);
      _toast('Optimize failed (see terminal).');
    } finally {
      if (mounted) setState(() => _busy = false);
    }
  }

  // ---------- Save/Load ----------
  Map<String, dynamic> _toJson() => {
    "version": 5,
    "poly": poly.map((p) => {"x": p.dx, "y": p.dy}).toList(),
    "speakers": speakers.map((p) => {"x": p.dx, "y": p.dy}).toList(),
    "seat": seat == null ? null : {"x": seat!.dx, "y": seat!.dy},
    "seatType": seatType.name,
    "allowedArea": allowedArea == null
        ? null
        : {
            "l": allowedArea!.left,
            "t": allowedArea!.top,
            "r": allowedArea!.right,
            "b": allowedArea!.bottom,
          },
    "settings": {
      "freqHz": freqHz,
      "optFreqs": optFreqs,
      "gain": gain,
      "wallRefl": wallRefl,
      "reflections": reflections,
      "gridStep": gridStep,
      "optSpeakers": optSpeakers,
      "optIters": optIters,
      "placementStyle": placementStyle.name,
      "robustHeatmapScale": robustHeatmapScale,
    },
  };

  void _fromJson(Map<String, dynamic> j) {
    final p = (j["poly"] as List)
        .cast<Map>()
        .map(
          (m) => Offset((m["x"] as num).toDouble(), (m["y"] as num).toDouble()),
        )
        .toList();

    final sp = (j["speakers"] as List)
        .cast<Map>()
        .map(
          (m) => Offset((m["x"] as num).toDouble(), (m["y"] as num).toDouble()),
        )
        .toList();

    Offset? s;
    final seatJ = j["seat"];
    if (seatJ != null) {
      s = Offset(
        (seatJ["x"] as num).toDouble(),
        (seatJ["y"] as num).toDouble(),
      );
    }

    Rect? aa;
    final aaJ = j["allowedArea"];
    if (aaJ != null) {
      aa = Rect.fromLTRB(
        (aaJ["l"] as num).toDouble(),
        (aaJ["t"] as num).toDouble(),
        (aaJ["r"] as num).toDouble(),
        (aaJ["b"] as num).toDouble(),
      );
    }

    final st = (j["seatType"] as String?) ?? "couch";
    final stEnum = (st == "chair") ? SeatType.chair : SeatType.couch;

    final sJ = (j["settings"] as Map?) ?? {};

    final psName = (sJ["placementStyle"] as String?) ?? placementStyle.name;
    final psEnum = PlacementStyle.values.firstWhere(
      (e) => e.name == psName,
      orElse: () => PlacementStyle.realistic,
    );

    setState(() {
      poly = p;
      speakers = sp;
      seat = s;
      seatType = stEnum;
      allowedArea = aa;

      freqHz = ((sJ["freqHz"] ?? freqHz) as num).toDouble();
      optFreqs = (sJ["optFreqs"] ?? optFreqs).toString();
      gain = ((sJ["gain"] ?? gain) as num).toDouble();
      wallRefl = ((sJ["wallRefl"] ?? wallRefl) as num).toDouble();
      reflections = (sJ["reflections"] ?? reflections) as int;
      gridStep = ((sJ["gridStep"] ?? gridStep) as num).toDouble();
      optSpeakers = (sJ["optSpeakers"] ?? optSpeakers) as int;
      optIters = (sJ["optIters"] ?? optIters) as int;
      placementStyle = psEnum;
      robustHeatmapScale =
          (sJ["robustHeatmapScale"] as bool?) ?? robustHeatmapScale;

      _heatImage?.dispose();
      _heatImage = null;
      _areaDragStartWorld = null;
    });

    _refreshViewSetState();
  }

  Future<void> _saveConfig() async {
    final jsonText = const JsonEncoder.withIndent("  ").convert(_toJson());
    final bytes = Uint8List.fromList(utf8.encode(jsonText));

    final dir = await getApplicationDocumentsDirectory();
    final fallbackPath = "${dir.path}/room_config.json";

    try {
      final picked = await FilePicker.platform.saveFile(
        dialogTitle: "Save configuration",
        fileName: "room_config.json",
        allowedExtensions: const ["json"],
        type: FileType.custom,
        bytes: bytes, // REQUIRED on Android/iOS
      );
      if (picked != null) {
        _toast("Saved.");
        return;
      }
    } catch (_) {
      // fall back below
    }

    await File(fallbackPath).writeAsString(jsonText);
    _toast("Saved to app storage.");
  }

  Future<void> _loadConfig() async {
    final picked = await FilePicker.platform.pickFiles(
      dialogTitle: "Load configuration",
      type: FileType.custom,
      allowedExtensions: const ["json"],
    );
    if (picked == null || picked.files.isEmpty) return;

    final path = picked.files.first.path;
    if (path == null) return;

    final txt = await File(path).readAsString();
    _fromJson(jsonDecode(txt) as Map<String, dynamic>);
    await _simulate();
    _toast("Loaded.");
  }

  // ---------- UI helpers ----------
  void _toast(String msg) {
    if (!mounted) return;
    ScaffoldMessenger.of(context).showSnackBar(
      SnackBar(content: Text(msg), duration: const Duration(seconds: 2)),
    );
  }

  Future<void> _openSettings() async {
    final res = await Navigator.of(context).push<SettingsResult>(
      MaterialPageRoute(
        builder: (_) => SettingsPage(
          initial: SettingsResult(
            freqHz: freqHz,
            optFreqs: optFreqs,
            gain: gain,
            wallRefl: wallRefl,
            reflections: reflections,
            gridStep: gridStep,
            optSpeakers: optSpeakers,
            optIters: optIters,
            seatType: seatType,
            occlusionEnabled: _occlusionEnabled,
            softDiffraction: _softDiffraction,
            placementStyle: placementStyle,
            robustHeatmapScale: robustHeatmapScale,
          ),
        ),
      ),
    );
    if (res == null) return;

    setState(() {
      freqHz = res.freqHz;
      optFreqs = res.optFreqs;
      gain = res.gain;
      wallRefl = res.wallRefl;
      reflections = res.reflections;
      gridStep = res.gridStep;
      optSpeakers = res.optSpeakers;
      optIters = res.optIters;
      seatType = res.seatType;
      _occlusionEnabled = res.occlusionEnabled;
      _softDiffraction = res.softDiffraction;
      placementStyle = res.placementStyle;
      robustHeatmapScale = res.robustHeatmapScale;
    });
  }

  // ---------- Tutorial ----------
  void _startTutorial() {
    setState(() {
      showTutorial = true;
      tutorialStep = 0;
    });
  }

  void _nextTutorialStep() {
    setState(() {
      tutorialStep += 1;
      if (tutorialStep > 5) showTutorial = false;
    });
  }

  void _skipTutorial() {
    setState(() {
      showTutorial = false;
      tutorialStep = 0;
    });
  }

  // ---------- Tap actions (called on TapUp only, so long-press can cancel it) ----------
  void _onTapAction(TapUpDetails d, Size size) {
    if (_busy) return;
    if (_ignoreNextTapUp) {
      _ignoreNextTapUp = false;
      return;
    }

    final world = _screenToWorld(d.localPosition, size, _viewWorld);

    if (mode != EditMode.room && !_pointInPoly(world, poly)) {
      _toast("Tap inside the room.");
      return;
    }

    setState(() {
      if (mode == EditMode.speakers) {
        speakers.add(world);
      } else if (mode == EditMode.seat) {
        seat = world;
      }
    });
  }

  // ---------- Build ----------
  @override
  Widget build(BuildContext context) {
    _refreshView();

    // Local coordinate frame for readouts & presets: origin at poly[0] (front-left),
    // X along poly[0]→poly[1] (front wall), Y along poly[0]→poly[3] (left wall into room).
    final _LocalFrame frame = _LocalFrame.fromPoly(poly);
    final PosRow? seatRow = (seat != null) ? frame.posRow('Seat', seat!) : null;
    final List<PosRow> speakerRows = [
      for (int i = 0; i < speakers.length; i++)
        frame.posRow('S${i + 1}', speakers[i]),
    ];

    return Scaffold(
      appBar: AppBar(
        title: const Text("Room Pressure Simulator"),
        actions: [
          IconButton(
            icon: const Icon(Icons.folder_open),
            onPressed: _busy ? null : _loadConfig,
            tooltip: "Load",
          ),
          IconButton(
            icon: const Icon(Icons.save),
            onPressed: _busy ? null : _saveConfig,
            tooltip: "Save",
          ),

          PopupMenuButton<StereoPreset>(
            tooltip: "Stereo presets",
            icon: const Icon(Icons.speaker),
            onSelected: (p) => _applyPreset(p),
            itemBuilder: (context) => const [
              PopupMenuItem(
                value: StereoPreset.typical,
                child: Text("Typical stereo"),
              ),
              PopupMenuItem(
                value: StereoPreset.wide,
                child: Text("Wide stereo"),
              ),
              PopupMenuItem(
                value: StereoPreset.tooWide,
                child: Text("Too wide (demo)"),
              ),
            ],
          ),
          IconButton(
            icon: const Icon(Icons.restart_alt),
            onPressed: _busy ? null : _resetToDefaults,
            tooltip: "Reset to defaults",
          ),
          IconButton(
            key: _kSettingsBtn,
            icon: const Icon(Icons.settings),
            onPressed: _busy ? null : _openSettings,
            tooltip: "Settings",
          ),
          IconButton(
            icon: const Icon(Icons.help_outline),
            tooltip: 'Help',
            onPressed: () {
              Navigator.of(context).push(
                MaterialPageRoute(
                  builder: (_) =>
                      IntroScreen(firstLaunch: false, onDone: () async {}),
                ),
              );
            },
          ),
          IconButton(
            icon: const Icon(Icons.school),
            onPressed: _startTutorial,
            tooltip: "Tutorial",
          ),
        ],
      ),
      body: Stack(
        children: [
          Column(
            children: [
              Expanded(
                child: LayoutBuilder(
                  builder: (context, constraints) {
                    final size = Size(
                      constraints.maxWidth,
                      constraints.maxHeight,
                    );

                    return Stack(
                      children: [
                        GestureDetector(
                          key: _kCanvas,

                          // Tap placement (speakers/seat)
                          onTapUp: (d) {
                            if (mode != EditMode.room &&
                                mode != EditMode.area) {
                              _onTapAction(d, size);
                            }
                          },

                          // Pan start: ROOM pick vertex, AREA start drag
                          onPanStart: (d) {
                            if (_busy) return;

                            if (mode == EditMode.room) {
                              final idx = _pickVertex(d.localPosition, size);
                              if (idx != null)
                                setState(() => draggingVertex = idx);
                              return;
                            }

                            if (mode == EditMode.area) {
                              final world = _screenToWorld(
                                d.localPosition,
                                size,
                                _viewWorld,
                              );
                              if (!_pointInPoly(world, poly)) return;
                              _ignoreNextTapUp = true;
                              _areaDragStartWorld = world;
                              setState(
                                () => allowedArea = Rect.fromLTRB(
                                  world.dx,
                                  world.dy,
                                  world.dx,
                                  world.dy,
                                ),
                              );
                            }
                          },

                          // Pan update: ROOM move vertex, AREA resize
                          onPanUpdate: (d) {
                            if (_busy) return;

                            if (mode == EditMode.room &&
                                draggingVertex != null) {
                              final w = _screenToWorld(
                                d.localPosition,
                                size,
                                _viewWorld,
                              );
                              setState(() => poly[draggingVertex!] = w);
                              return;
                            }

                            if (mode == EditMode.area &&
                                _areaDragStartWorld != null) {
                              final world = _screenToWorld(
                                d.localPosition,
                                size,
                                _viewWorld,
                              );
                              final a = _areaDragStartWorld!;
                              setState(() {
                                allowedArea = Rect.fromLTRB(
                                  math.min(a.dx, world.dx),
                                  math.min(a.dy, world.dy),
                                  math.max(a.dx, world.dx),
                                  math.max(a.dy, world.dy),
                                );
                              });
                            }
                          },

                          onPanEnd: (_) {
                            if (mode == EditMode.room) {
                              setState(() => draggingVertex = null);
                              _refreshViewSetState();
                              return;
                            }
                            if (mode == EditMode.area) {
                              _areaDragStartWorld = null;
                            }
                          },

                          // ROOM: add corner
                          onDoubleTapDown: (d) {
                            if (_busy) return;
                            if (mode != EditMode.room) return;
                            final edge = _pickEdgeForInsert(
                              d.localPosition,
                              size,
                            );
                            if (edge == null) return;
                            final w = _screenToWorld(
                              d.localPosition,
                              size,
                              _viewWorld,
                            );
                            _insertCornerAtEdge(edge, w);
                            _refreshViewSetState();
                          },

                          // Long press: ROOM delete, SPEAKERS remove, AREA clear
                          onLongPressStart: (d) {
                            if (_busy) return;

                            if (mode == EditMode.room) {
                              final idx = _pickVertex(d.localPosition, size);
                              if (idx != null) {
                                _deleteCorner(idx);
                                _refreshViewSetState();
                                _ignoreNextTapUp = true;
                              }
                              return;
                            }

                            if (mode == EditMode.speakers &&
                                speakers.isNotEmpty) {
                              final world = _screenToWorld(
                                d.localPosition,
                                size,
                                _viewWorld,
                              );
                              int best = 0;
                              double bestD = double.infinity;
                              for (int i = 0; i < speakers.length; i++) {
                                final d2 =
                                    (speakers[i] - world).distanceSquared;
                                if (d2 < bestD) {
                                  bestD = d2;
                                  best = i;
                                }
                              }
                              setState(() => speakers.removeAt(best));
                              _ignoreNextTapUp = true;
                              return;
                            }

                            if (mode == EditMode.area && allowedArea != null) {
                              final world = _screenToWorld(
                                d.localPosition,
                                size,
                                _viewWorld,
                              );
                              if (allowedArea!.contains(world)) {
                                setState(() {
                                  allowedArea = null;
                                  _areaDragStartWorld = null;
                                });
                                _ignoreNextTapUp = true;
                                _toast("Area cleared.");
                              }
                            }
                          },

                          child: CustomPaint(
                            painter: RoomPainter(
                              poly: poly,
                              speakers: speakers,
                              seat: seat,
                              seatType: seatType,
                              allowedArea: allowedArea,
                              mode: mode,
                              viewWorld: _viewWorld,
                              heatWorld: _heatWorld,
                              heatImage: _heatImage,
                              maxSpl: _lastMaxSpl,
                            ),
                            child: const SizedBox.expand(),
                          ),
                        ),
                      ],
                    );
                  },
                ),
              ),

              // Legend sits BELOW the map (never blocks gestures / never covers the heatmap).
              if (_heatImage != null)
                Padding(
                  padding: const EdgeInsets.fromLTRB(12, 10, 12, 0),
                  child: RoomHeatmapLegend(
                    hasHeat: _heatImage != null,
                    seatRelDb: _lastSeatRelDb,
                    robustScaleUsed: robustHeatmapScale,
                    seatRow: seatRow,
                    speakerRows: speakerRows,
                    showDetails: _showLegendDetails,
                    onShowDetailsChanged: (v) {
                      setState(() => _showLegendDetails = v);
                    },
                  ),
                ),
              const SizedBox(height: 10),
              SafeArea(
                top: false,
                child: Padding(
                  padding: EdgeInsets.fromLTRB(
                    12,
                    8,
                    12,
                    10 + MediaQuery.of(context).padding.bottom,
                  ),
                  child: Column(
                    children: [
                      Row(
                        key: _kModeRow,
                        children: [
                          _modeChip(label: "Room", value: EditMode.room),
                          const SizedBox(width: 8),
                          _modeChip(
                            label: "Speakers",
                            value: EditMode.speakers,
                          ),
                          const SizedBox(width: 8),
                          _modeChip(label: "Seat", value: EditMode.seat),
                          const SizedBox(width: 8),
                          _modeChip(label: "Area", value: EditMode.area),
                        ],
                      ),
                      const SizedBox(height: 10),
                      Row(
                        children: [
                          Expanded(
                            child: FilledButton.icon(
                              key: _kSimBtn,
                              onPressed: _busy ? null : _simulate,
                              icon: const Icon(Icons.play_arrow),
                              label: Text(_busy ? "Working..." : "Simulate"),
                            ),
                          ),
                          const SizedBox(width: 10),
                          Expanded(
                            child: FilledButton.icon(
                              key: _kOptBtn,
                              onPressed: _busy ? null : _optimize,
                              icon: const Icon(Icons.auto_fix_high),
                              label: const Text("Optimize"),
                            ),
                          ),
                        ],
                      ),
                      const SizedBox(height: 8),
                      Row(
                        children: [
                          Expanded(
                            child: Text(
                              "Mode: ${mode.name.toUpperCase()}  |  Speakers: ${speakers.length}  |  Seat: ${seat == null ? 'no' : 'yes'}",
                              style: Theme.of(context).textTheme.bodySmall,
                            ),
                          ),
                          IconButton(
                            onPressed: () {
                              setState(() {
                                _heatImage?.dispose();
                                _heatImage = null;
                                _lastMaxSpl = null;
                                _lastSeatRelDb = null;
                              });
                            },
                            tooltip: "Clear heatmap",
                            icon: const Icon(Icons.layers_clear),
                          ),
                        ],
                      ),
                    ],
                  ),
                ),
              ),
            ],
          ),

          if (showTutorial)
            TutorialOverlay(
              step: tutorialStep,
              onNext: _nextTutorialStep,
              onSkip: _skipTutorial,
              anchors: TutorialAnchors(
                modeRow: _kModeRow,
                canvas: _kCanvas,
                simulate: _kSimBtn,
                optimize: _kOptBtn,
                settings: _kSettingsBtn,
              ),
            ),

          // Loading overlay
          if (_busy)
            Positioned.fill(
              child: AbsorbPointer(
                child: Container(
                  color: Colors.black.withValues(alpha: 0.35),
                  child: Center(
                    child: Card(
                      elevation: 6,
                      child: Padding(
                        padding: const EdgeInsets.symmetric(
                          horizontal: 24,
                          vertical: 18,
                        ),
                        child: Column(
                          mainAxisSize: MainAxisSize.min,
                          children: const [
                            CircularProgressIndicator(),
                            SizedBox(height: 14),
                            Text("Working…", style: TextStyle(fontSize: 14)),
                          ],
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            ),
        ],
      ),
    );
  }

  Widget _modeChip({required String label, required EditMode value}) {
    return ChoiceChip(
      label: Text(label),
      selected: mode == value,
      onSelected: (_) {
        setState(() {
          mode = value;
          draggingVertex = null;
          _ignoreNextTapUp = false;
          _areaDragStartWorld = null;
        });
      },
    );
  }
}

const floatNaN = double.nan;

// ------------------------------------------------------------
// Painter
// ------------------------------------------------------------
class RoomPainter extends CustomPainter {
  final List<Offset> poly;
  final List<Offset> speakers;
  final Offset? seat;
  final SeatType seatType;
  final Rect? allowedArea;
  final EditMode mode;

  final Rect viewWorld;
  final Rect heatWorld;
  final ui.Image? heatImage;
  final double? maxSpl;

  RoomPainter({
    required this.poly,
    required this.speakers,
    required this.seat,
    required this.seatType,
    required this.allowedArea,
    required this.mode,
    required this.viewWorld,
    required this.heatWorld,
    required this.heatImage,
    required this.maxSpl,
  });

  Offset w2s(Offset w, Size size) {
    final sx = (w.dx - viewWorld.left) / viewWorld.width * size.width;
    final sy = (w.dy - viewWorld.top) / viewWorld.height * size.height;
    return Offset(sx, sy);
  }

  @override
  void paint(Canvas canvas, Size size) {
    // Heatmap
    if (heatImage != null) {
      final src = Rect.fromLTWH(
        0,
        0,
        heatImage!.width.toDouble(),
        heatImage!.height.toDouble(),
      );

      final topLeft = w2s(heatWorld.topLeft, size);
      final bottomRight = w2s(heatWorld.bottomRight, size);

      final dst = Rect.fromLTRB(
        math.min(topLeft.dx, bottomRight.dx),
        math.min(topLeft.dy, bottomRight.dy),
        math.max(topLeft.dx, bottomRight.dx),
        math.max(topLeft.dy, bottomRight.dy),
      );

      canvas.drawImageRect(
        heatImage!,
        src,
        dst,
        Paint()..filterQuality = FilterQuality.none,
      );
    }

    // Room outline
    final roomPath = Path();
    if (poly.isNotEmpty) {
      final p0 = w2s(poly.first, size);
      roomPath.moveTo(p0.dx, p0.dy);
      for (int i = 1; i < poly.length; i++) {
        final p = w2s(poly[i], size);
        roomPath.lineTo(p.dx, p.dy);
      }
      roomPath.close();

      canvas.drawPath(
        roomPath,
        Paint()
          ..style = PaintingStyle.stroke
          ..strokeWidth = 3
          ..color = Colors.blueGrey,
      );
    }

    // Allowed area
    if (allowedArea != null) {
      final r = allowedArea!;
      final a = w2s(Offset(r.left, r.top), size);
      final b = w2s(Offset(r.right, r.bottom), size);
      final rect = Rect.fromLTRB(
        math.min(a.dx, b.dx),
        math.min(a.dy, b.dy),
        math.max(a.dx, b.dx),
        math.max(a.dy, b.dy),
      );
      _drawDashedRect(canvas, rect, Colors.black, 2);
    }

    // Wall labels + angles
    _drawWallLabels(canvas, size);

    // Corners visible (room mode)
    if (mode == EditMode.room) {
      final paintV = Paint()..color = Colors.black;
      for (final v in poly) {
        final s = w2s(v, size);
        canvas.drawCircle(s, 7, paintV);
      }
    }

    // Speakers
    final spPaint = Paint()..color = Colors.red;
    for (int i = 0; i < speakers.length; i++) {
      final s = w2s(speakers[i], size);
      _drawSpeakerIcon(canvas, s, spPaint);
      final tp = TextPainter(
        text: TextSpan(
          text: "S${i + 1}",
          style: const TextStyle(
            color: Colors.black,
            fontSize: 12,
            fontWeight: FontWeight.w700,
          ),
        ),
        textDirection: TextDirection.ltr,
      )..layout();
      tp.paint(canvas, s + const Offset(10, -10));
    }

    // Seat (simple shape)
    if (seat != null) {
      final s = w2s(seat!, size);
      final seatPaint = Paint()..color = Colors.deepPurple;
      if (seatType == SeatType.couch) {
        final rect = RRect.fromRectAndRadius(
          Rect.fromCenter(center: s, width: 44, height: 22),
          const Radius.circular(8),
        );
        canvas.drawRRect(rect, seatPaint);
      } else {
        final rect = RRect.fromRectAndRadius(
          Rect.fromCenter(center: s, width: 24, height: 24),
          const Radius.circular(6),
        );
        canvas.drawRRect(rect, seatPaint);
      }

      final lab = TextPainter(
        text: const TextSpan(
          text: "Seat",
          style: TextStyle(
            color: Colors.black,
            fontSize: 12,
            fontWeight: FontWeight.w600,
          ),
        ),
        textDirection: TextDirection.ltr,
      )..layout();
      lab.paint(canvas, s + const Offset(-14, 18));
    }

    // Hint
    final hint = _hintText();
    final ht = TextPainter(
      text: TextSpan(
        text: hint,
        style: TextStyle(
          color: Colors.black.withValues(alpha: 0.7),
          fontSize: 12,
        ),
      ),
      textDirection: TextDirection.ltr,
    )..layout(maxWidth: size.width - 16);
    ht.paint(canvas, const Offset(8, 8));
  }

  void _drawSpeakerIcon(Canvas canvas, Offset c, Paint basePaint) {
    // Simple "speaker cabinet" symbol (rounded rectangle + two drivers).
    const double w = 14;
    const double h = 20;
    final rect = Rect.fromCenter(center: c, width: w, height: h);
    final rrect = RRect.fromRectAndRadius(rect, const Radius.circular(3));

    final fill = Paint()
      ..style = PaintingStyle.fill
      ..color = Colors.white;
    final stroke = Paint()
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2
      ..color = basePaint.color;

    canvas.drawRRect(rrect, fill);
    canvas.drawRRect(rrect, stroke);

    final tweeter = Paint()..color = basePaint.color.withValues(alpha: 0.85);
    final woofer = Paint()..color = basePaint.color.withValues(alpha: 0.85);

    canvas.drawCircle(c.translate(0, -5), 2.2, tweeter);
    canvas.drawCircle(c.translate(0, 4.5), 3.6, woofer);
  }

  String _hintText() {
    switch (mode) {
      case EditMode.room:
        return "Room: drag corners. Double-tap edge to add. Long-press corner to delete.";
      case EditMode.speakers:
        return "Speakers: tap to add. Long-press near a speaker to remove.";
      case EditMode.seat:
        return "Seat: tap to place the seat.";
      case EditMode.area:
        return "Area: drag to draw allowed speaker area. Long-press inside it to clear.";
    }
  }

  void _drawWallLabels(Canvas canvas, Size size) {
    if (poly.length < 3) return;

    // Use centroid to decide inside/outside directions.
    Offset centroid = Offset.zero;
    for (final p in poly) {
      centroid += p;
    }
    centroid = centroid / poly.length.toDouble();

    final diag = math.sqrt(
      viewWorld.width * viewWorld.width + viewWorld.height * viewWorld.height,
    );
    final edgeOffset = (0.05 * diag).clamp(
      0.14,
      0.28,
    ); // pushed outside for placement clarity
    final cornerOffset = (0.06 * diag).clamp(
      0.18,
      0.36,
    ); // pushed outside for placement clarity

    // -------- Length labels aligned with walls --------
    for (int i = 0; i < poly.length; i++) {
      final a = poly[i];
      final b = poly[(i + 1) % poly.length];
      final e = b - a;
      final len = e.distance;
      if (len < 1e-9) continue;

      final mid = (a + b) / 2.0;

      // Outward normal: pick the normal pointing away from centroid.
      final nRaw = Offset(e.dy, -e.dx) / len; // one of the normals
      final toC = centroid - mid;
      final nOut = (toC.dx * nRaw.dx + toC.dy * nRaw.dy) > 0 ? -nRaw : nRaw;

      final labelPos = mid + nOut * edgeOffset;

      final txt = "${len.toStringAsFixed(2)} m";
      final tp = TextPainter(
        text: TextSpan(
          text: txt,
          style: const TextStyle(
            fontSize: 12,
            fontWeight: FontWeight.w700,
            color: Colors.black,
          ),
        ),
        textDirection: TextDirection.ltr,
      )..layout();

      final ang = math.atan2(e.dy, e.dx);

      // Background pill in local coords (rotated with the edge).
      final padX = 10.0;
      final padY = 6.0;
      final rect = Rect.fromCenter(
        center: Offset.zero,
        width: tp.width + padX * 2,
        height: tp.height + padY * 2,
      );
      final rrect = RRect.fromRectAndRadius(rect, const Radius.circular(10));

      canvas.save();
      final s = w2s(labelPos, size);
      canvas.translate(s.dx, s.dy);
      canvas.rotate(ang);

      canvas.drawRRect(
        rrect,
        Paint()
          ..color = Colors.white.withValues(alpha: 0.85)
          ..style = PaintingStyle.fill,
      );
      canvas.drawRRect(
        rrect,
        Paint()
          ..style = PaintingStyle.stroke
          ..strokeWidth = 1
          ..color = Colors.black.withValues(alpha: 0.35),
      );

      tp.paint(canvas, Offset(-tp.width / 2, -tp.height / 2));
      canvas.restore();
    }

    // -------- Corner angles anchored at corners --------
    for (int i = 0; i < poly.length; i++) {
      final prev = poly[(i - 1 + poly.length) % poly.length];
      final cur = poly[i];
      final next = poly[(i + 1) % poly.length];

      final v1 = prev - cur;
      final v2 = next - cur;
      final n1 = v1.distance;
      final n2 = v2.distance;
      if (n1 < 1e-9 || n2 < 1e-9) continue;

      final dot = (v1.dx * v2.dx + v1.dy * v2.dy) / (n1 * n2);
      final angDeg = (math.acos(dot.clamp(-1.0, 1.0)) * 180 / math.pi);

      final u1 = v1 / n1;
      final u2 = v2 / n2;
      var bis = u1 + u2;
      if (bis.distance < 1e-9) {
        bis = Offset(-u1.dy, u1.dx); // perpendicular fallback
      } else {
        bis = bis / bis.distance;
      }

      // Ensure bisector points *outward* (away from centroid) so labels never cover the interior.
      final toC = centroid - cur; // points inward
      if ((toC.dx * bis.dx + toC.dy * bis.dy) > 0) bis = -bis;

      final p = cur + bis * cornerOffset;

      final txt = "${angDeg.round()}°";
      final tp = TextPainter(
        text: TextSpan(
          text: txt,
          style: const TextStyle(
            fontSize: 12,
            fontWeight: FontWeight.w800,
            color: Colors.black,
          ),
        ),
        textDirection: TextDirection.ltr,
      )..layout();

      final padX = 10.0;
      final padY = 6.0;
      final rect = Rect.fromCenter(
        center: Offset.zero,
        width: tp.width + padX * 2,
        height: tp.height + padY * 2,
      );
      final rrect = RRect.fromRectAndRadius(rect, const Radius.circular(10));

      canvas.save();
      final s = w2s(p, size);
      canvas.translate(s.dx, s.dy);

      canvas.drawRRect(
        rrect,
        Paint()
          ..color = Colors.white.withValues(alpha: 0.85)
          ..style = PaintingStyle.fill,
      );
      canvas.drawRRect(
        rrect,
        Paint()
          ..style = PaintingStyle.stroke
          ..strokeWidth = 1
          ..color = Colors.black.withValues(alpha: 0.35),
      );
      tp.paint(canvas, Offset(-tp.width / 2, -tp.height / 2));
      canvas.restore();
    }
  }

  void _drawDashedRect(Canvas canvas, Rect rect, Color color, double w) {
    final p = Paint()
      ..color = color
      ..strokeWidth = w
      ..style = PaintingStyle.stroke;
    const dash = 8.0;
    const gap = 5.0;

    void drawDashedLine(Offset a, Offset b) {
      final dx = b.dx - a.dx;
      final dy = b.dy - a.dy;
      final len = math.sqrt(dx * dx + dy * dy);
      if (len < 1e-6) return;
      final ux = dx / len;
      final uy = dy / len;
      double t = 0;
      while (t < len) {
        final t2 = math.min(t + dash, len);
        canvas.drawLine(
          Offset(a.dx + ux * t, a.dy + uy * t),
          Offset(a.dx + ux * t2, a.dy + uy * t2),
          p,
        );
        t += dash + gap;
      }
    }

    drawDashedLine(rect.topLeft, rect.topRight);
    drawDashedLine(rect.topRight, rect.bottomRight);
    drawDashedLine(rect.bottomRight, rect.bottomLeft);
    drawDashedLine(rect.bottomLeft, rect.topLeft);
  }

  @override
  bool shouldRepaint(covariant RoomPainter old) {
    return old.poly != poly ||
        old.speakers != speakers ||
        old.seat != seat ||
        old.seatType != seatType ||
        old.allowedArea != allowedArea ||
        old.mode != mode ||
        old.viewWorld != viewWorld ||
        old.heatImage != heatImage ||
        old.heatWorld != heatWorld ||
        old.maxSpl != maxSpl;
  }
}

// ------------------------------------------------------------
// Tutorial Overlay
// ------------------------------------------------------------
class TutorialAnchors {
  final GlobalKey modeRow;
  final GlobalKey canvas;
  final GlobalKey simulate;
  final GlobalKey optimize;
  final GlobalKey settings;
  const TutorialAnchors({
    required this.modeRow,
    required this.canvas,
    required this.simulate,
    required this.optimize,
    required this.settings,
  });
}

class TutorialOverlay extends StatelessWidget {
  final int step;
  final VoidCallback onNext;
  final VoidCallback onSkip;
  final TutorialAnchors anchors;

  const TutorialOverlay({
    super.key,
    required this.step,
    required this.onNext,
    required this.onSkip,
    required this.anchors,
  });

  Rect? _rectOf(GlobalKey key) {
    final ctx = key.currentContext;
    if (ctx == null) return null;
    final box = ctx.findRenderObject() as RenderBox?;
    if (box == null) return null;
    final pos = box.localToGlobal(Offset.zero);
    return pos & box.size;
  }

  @override
  Widget build(BuildContext context) {
    final rMode = _rectOf(anchors.modeRow);
    final rCanvas = _rectOf(anchors.canvas);
    final rSim = _rectOf(anchors.simulate);
    final rOpt = _rectOf(anchors.optimize);

    Rect? focus;
    String title = "";
    String body = "";
    bool showNext = true;

    switch (step) {
      case 0:
        focus = rMode;
        title = "Step 1: Choose a mode";
        body =
            "Use these buttons to switch between editing the room, placing speakers, placing the seat, and drawing an allowed area.";
        break;
      case 1:
        focus = rCanvas;
        title = "Step 2: Edit the room";
        body =
            "In ROOM mode:\n• Drag corners\n• Double-tap an edge to add a corner\n• Long-press a corner to delete it";
        break;
      case 2:
        focus = rMode;
        title = "Step 3: Place speakers";
        body =
            "Switch to SPEAKERS mode:\n• Tap to add speakers\n• Long-press near a speaker to remove it";
        break;
      case 3:
        focus = rMode;
        title = "Step 4: Place your seat";
        body =
            "Switch to SEAT mode and tap where you sit. Optimization uses this point.";
        break;
      case 4:
        focus = rMode;
        title = "Step 5: Allowed area (optional)";
        body =
            "Switch to AREA mode and DRAG to draw where speakers are allowed.\nLong-press inside the area to clear it.";
        break;
      case 5:
        focus = rSim ?? rOpt;
        title = "Step 6: Simulate and optimize";
        body =
            "Tap SIMULATE to draw the heatmap for the current frequency.\nTap OPTIMIZE to find speaker positions that are smoothest at the seat across Opt freqs (in Settings).";
        showNext = false;
        break;
      default:
        showNext = false;
    }

    final media = MediaQuery.of(context);
    final full = Rect.fromLTWH(0, 0, media.size.width, media.size.height);

    return Material(
      color: Colors.black.withValues(alpha: 0.55),
      child: Stack(
        children: [
          if (focus != null)
            Positioned.fromRect(
              rect: focus.inflate(6).intersect(full),
              child: IgnorePointer(
                child: Container(
                  decoration: BoxDecoration(
                    borderRadius: BorderRadius.circular(14),
                    border: Border.all(color: Colors.white, width: 3),
                  ),
                ),
              ),
            ),
          Align(
            alignment: Alignment.bottomCenter,
            child: Padding(
              padding: const EdgeInsets.fromLTRB(12, 12, 12, 18),
              child: Card(
                child: Padding(
                  padding: const EdgeInsets.fromLTRB(16, 14, 16, 12),
                  child: Column(
                    mainAxisSize: MainAxisSize.min,
                    children: [
                      Row(
                        children: [
                          Expanded(
                            child: Text(
                              title,
                              style: Theme.of(context).textTheme.titleMedium,
                            ),
                          ),
                          IconButton(
                            onPressed: onSkip,
                            icon: const Icon(Icons.close),
                            tooltip: "Close tutorial",
                          ),
                        ],
                      ),
                      const SizedBox(height: 6),
                      Align(alignment: Alignment.centerLeft, child: Text(body)),
                      const SizedBox(height: 12),
                      Row(
                        children: [
                          TextButton(
                            onPressed: onSkip,
                            child: const Text("Skip"),
                          ),
                          const Spacer(),
                          FilledButton(
                            onPressed: showNext ? onNext : onSkip,
                            child: Text(showNext ? "Next" : "Done"),
                          ),
                        ],
                      ),
                    ],
                  ),
                ),
              ),
            ),
          ),
        ],
      ),
    );
  }
}

// ------------------------------------------------------------
// Settings Page
// ------------------------------------------------------------
class SettingsResult {
  final double freqHz;
  final String optFreqs;
  final double gain;
  final double wallRefl;
  final int reflections;
  final double gridStep;
  final int optSpeakers;
  final int optIters;
  final SeatType seatType;
  final bool occlusionEnabled;
  final bool softDiffraction;
  final PlacementStyle placementStyle;
  final bool robustHeatmapScale;

  const SettingsResult({
    required this.freqHz,
    required this.optFreqs,
    required this.gain,
    required this.wallRefl,
    required this.reflections,
    required this.gridStep,
    required this.optSpeakers,
    required this.optIters,
    required this.seatType,
    required this.occlusionEnabled,
    required this.softDiffraction,
    required this.placementStyle,
    required this.robustHeatmapScale,
  });

  SettingsResult copyWith({
    double? freqHz,
    String? optFreqs,
    double? gain,
    double? wallRefl,
    int? reflections,
    double? gridStep,
    int? optSpeakers,
    int? optIters,
    SeatType? seatType,
    bool? occlusionEnabled,
    bool? softDiffraction,
    PlacementStyle? placementStyle,
    bool? robustHeatmapScale,
  }) {
    return SettingsResult(
      freqHz: freqHz ?? this.freqHz,
      optFreqs: optFreqs ?? this.optFreqs,
      gain: gain ?? this.gain,
      wallRefl: wallRefl ?? this.wallRefl,
      reflections: reflections ?? this.reflections,
      gridStep: gridStep ?? this.gridStep,
      optSpeakers: optSpeakers ?? this.optSpeakers,
      optIters: optIters ?? this.optIters,
      seatType: seatType ?? this.seatType,
      occlusionEnabled: occlusionEnabled ?? this.occlusionEnabled,
      softDiffraction: softDiffraction ?? this.softDiffraction,
      placementStyle: placementStyle ?? this.placementStyle,
      robustHeatmapScale: robustHeatmapScale ?? this.robustHeatmapScale,
    );
  }
}

class _LocalFrame {
  _LocalFrame(this.origin, this.xAxis, this.yAxis);

  final Offset origin;
  final Offset xAxis; // unit
  final Offset yAxis; // unit

  static _LocalFrame fromPoly(List<Offset> poly) {
    // Fallback to axis-aligned 1m frame if geometry is degenerate.
    if (poly.length < 4) {
      return _LocalFrame(
        const Offset(0, 0),
        const Offset(1, 0),
        const Offset(0, 1),
      );
    }
    final o = poly[0];
    final vx = poly[1] - poly[0];
    final vy = poly[3] - poly[0];

    Offset norm(Offset v) {
      final l = v.distance;
      if (l < 1e-9) return const Offset(0, 0);
      return Offset(v.dx / l, v.dy / l);
    }

    final x = norm(vx);
    final y = norm(vy);

    return _LocalFrame(o, x, y);
  }

  Offset toLocal(Offset world) {
    final d = world - origin;
    final x = d.dx * xAxis.dx + d.dy * xAxis.dy;
    final y = d.dx * yAxis.dx + d.dy * yAxis.dy;
    return Offset(x, y);
  }

  PosRow posRow(String label, Offset world) {
    final p = toLocal(world);
    return PosRow(label: label, x: p.dx, y: p.dy);
  }
}

class PosRow {
  const PosRow({required this.label, required this.x, required this.y});
  final String label;
  final double x;
  final double y;
}

class RoomHeatmapLegend extends StatelessWidget {
  const RoomHeatmapLegend({
    super.key,
    required this.hasHeat,
    required this.seatRelDb,
    this.seatRow,
    this.speakerRows = const [],
    this.robustScaleUsed = false,
    this.showDetails = true,
    this.onShowDetailsChanged,
  });

  final bool hasHeat;
  final double? seatRelDb;
  final PosRow? seatRow;
  final List<PosRow> speakerRows;

  /// True when the heatmap is normalized to a robust percentile reference (instead of strict max).
  final bool robustScaleUsed;

  /// When false, the legend collapses to the most compact form (gradient + ticks).
  final bool showDetails;

  /// If provided, shows an expand/collapse control.
  final ValueChanged<bool>? onShowDetailsChanged;

  // Must match the heatmap normalization: (spl - zRef).clamp(-30, 0)
  static const LinearGradient _gradient = LinearGradient(
    begin: Alignment.centerLeft,
    end: Alignment.centerRight,
    colors: [
      Color.fromARGB(255, 0, 176, 0),
      Color.fromARGB(255, 255, 255, 0),
      Color.fromARGB(255, 255, 0, 0),
    ],
    stops: [0.0, 0.5, 1.0],
  );

  @override
  Widget build(BuildContext context) {
    final bool valid = hasHeat;

    final TextStyle subStyle = TextStyle(
      fontSize: 11,
      color: Colors.black.withValues(alpha: 0.70),
      fontWeight: FontWeight.w600,
    );

    String fmt1(double v) => v.toStringAsFixed(1);

    Widget gradientBar({required double height}) {
      return ClipRRect(
        borderRadius: BorderRadius.circular(7),
        child: SizedBox(
          height: height,
          width: double.infinity,
          child: const DecoratedBox(
            decoration: BoxDecoration(gradient: _gradient),
          ),
        ),
      );
    }

    Widget ticksRow({required bool compact}) {
      return Row(
        mainAxisAlignment: MainAxisAlignment.spaceBetween,
        children: [
          _LegendTick(label: '-30 dB', compact: compact),
          _LegendTick(label: '-15 dB', compact: compact),
          _LegendTick(label: '0 dB', compact: compact),
        ],
      );
    }

    Widget expandCollapseButton({required bool expand}) {
      if (onShowDetailsChanged == null) return const SizedBox.shrink();
      return IconButton(
        tooltip: expand ? 'Show details' : 'Hide details',
        padding: EdgeInsets.zero,
        constraints: const BoxConstraints(minWidth: 28, minHeight: 28),
        onPressed: () => onShowDetailsChanged!.call(expand),
        icon: Icon(expand ? Icons.expand_more : Icons.expand_less, size: 18),
      );
    }

    // Collapsed: as compact as possible while retaining function/readability.
    Widget collapsed() {
      return Padding(
        padding: const EdgeInsets.fromLTRB(10, 6, 6, 6),
        child: Column(
          mainAxisSize: MainAxisSize.min,
          children: [
            Row(
              crossAxisAlignment: CrossAxisAlignment.center,
              children: [
                Expanded(child: gradientBar(height: 12)),
                const SizedBox(width: 6),
                expandCollapseButton(expand: true),
              ],
            ),
            const SizedBox(height: 4),
            ticksRow(compact: true),
          ],
        ),
      );
    }

    // Expanded: full legend + (optional) seat readout and positions.
    Widget expanded() {
      final rowsAvailable = (seatRow != null) || speakerRows.isNotEmpty;

      return Padding(
        padding: const EdgeInsets.fromLTRB(12, 10, 12, 10),
        child: ConstrainedBox(
          constraints: const BoxConstraints(minWidth: 200, maxWidth: 320),
          child: Column(
            mainAxisSize: MainAxisSize.min,
            crossAxisAlignment: CrossAxisAlignment.start,
            children: [
              Row(
                children: [
                  const Expanded(
                    child: Text(
                      'Heatmap',
                      style: TextStyle(
                        fontSize: 12,
                        fontWeight: FontWeight.w800,
                      ),
                    ),
                  ),
                  if (onShowDetailsChanged != null)
                    expandCollapseButton(expand: false),
                  if (rowsAvailable) ...[
                    const SizedBox(width: 6),
                    IconButton(
                      tooltip: 'Copy positions',
                      padding: EdgeInsets.zero,
                      constraints: const BoxConstraints(
                        minWidth: 28,
                        minHeight: 28,
                      ),
                      onPressed: () {
                        final buf = StringBuffer();
                        if (seatRow != null) buf.writeln(_formatRow(seatRow!));
                        for (final r in speakerRows) {
                          buf.writeln(_formatRow(r));
                        }
                        if (buf.isNotEmpty) {
                          Clipboard.setData(
                            ClipboardData(text: buf.toString().trim()),
                          );
                          ScaffoldMessenger.of(context).showSnackBar(
                            const SnackBar(
                              content: Text('Copied positions to clipboard'),
                            ),
                          );
                        }
                      },
                      icon: const Icon(Icons.copy, size: 18),
                    ),
                  ],
                ],
              ),
              const SizedBox(height: 6),
              Text(
                valid
                    ? (robustScaleUsed
                          ? 'Relative scale (0 dB = 98th percentile)'
                          : 'Relative scale (0 dB = max point)')
                    : 'Run Simulate to show heatmap',
                style: subStyle,
              ),
              const SizedBox(height: 10),
              gradientBar(height: 14),
              const SizedBox(height: 8),
              ticksRow(compact: false),
              if (showDetails) ...[
                const SizedBox(height: 10),
                if (valid && seatRelDb != null && seatRelDb!.isFinite)
                  Text('Seat: ${fmt1(seatRelDb!)} dB', style: subStyle),
                if (rowsAvailable) ...[
                  const SizedBox(height: 8),
                  _PosTable(seatRow: seatRow, speakerRows: speakerRows),
                ] else ...[
                  const SizedBox(height: 8),
                  Text(
                    'Place seat and speakers to see coordinates',
                    style: subStyle,
                  ),
                ],
              ],
            ],
          ),
        ),
      );
    }

    return DecoratedBox(
      decoration: BoxDecoration(
        color: Colors.white.withValues(alpha: 0.92),
        borderRadius: BorderRadius.circular(12),
        border: Border.all(color: Colors.black.withValues(alpha: 0.12)),
        boxShadow: [
          BoxShadow(
            blurRadius: 8,
            offset: const Offset(0, 2),
            color: Colors.black.withValues(alpha: 0.12),
          ),
        ],
      ),
      child: showDetails ? expanded() : collapsed(),
    );
  }
}

class _LegendTick extends StatelessWidget {
  const _LegendTick({required this.label, this.compact = false});

  final String label;
  final bool compact;

  @override
  Widget build(BuildContext context) {
    return Text(
      label,
      style: TextStyle(
        fontSize: compact ? 10.5 : 11,
        height: 1.0,
        color: Colors.black.withValues(alpha: 0.80),
        fontWeight: FontWeight.w700,
      ),
    );
  }
}

String _formatRow(PosRow r) {
  // Tabular-friendly row: "Seat  x=1.23m  y=2.34m"
  String f(double v) => v.toStringAsFixed(2).padLeft(6);
  return '${r.label.padRight(4)}  x=${f(r.x)} m  y=${f(r.y)} m';
}

class _PosTable extends StatelessWidget {
  const _PosTable({this.seatRow, required this.speakerRows});

  final PosRow? seatRow;
  final List<PosRow> speakerRows;

  @override
  Widget build(BuildContext context) {
    final rows = <PosRow>[if (seatRow != null) seatRow!, ...speakerRows];

    // Use tabular figures to align numbers nicely.
    final style = const TextStyle(
      fontSize: 11,
      fontWeight: FontWeight.w700,
      fontFeatures: [FontFeature.tabularFigures()],
    );

    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        for (final r in rows) ...[
          Text(_formatRow(r), style: style),
          const SizedBox(height: 2),
        ],
      ],
    );
  }
}

class SettingsPage extends StatefulWidget {
  final SettingsResult initial;
  const SettingsPage({super.key, required this.initial});

  @override
  State<SettingsPage> createState() => _SettingsPageState();
}

class _SettingsPageState extends State<SettingsPage> {
  late SettingsResult s;
  late TextEditingController _optFreqCtrl;

  @override
  void initState() {
    super.initState();
    s = widget.initial;
    _optFreqCtrl = TextEditingController(text: s.optFreqs);
  }

  @override
  void dispose() {
    _optFreqCtrl.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text("Settings"),
        actions: [
          TextButton(
            onPressed: () => Navigator.of(
              context,
            ).pop(s.copyWith(optFreqs: _optFreqCtrl.text)),
            child: const Text("Done"),
          ),
        ],
      ),
      body: ListView(
        padding: const EdgeInsets.fromLTRB(16, 12, 16, 24),
        children: [
          _sectionTitle(context, "Simulation"),
          _slider(
            "Single frequency (Hz)",
            s.freqHz,
            20,
            120,
            (v) => setState(() => s = s.copyWith(freqHz: v)),
            (v) => v.toStringAsFixed(0),
          ),
          _slider(
            "Grid step (m) (smaller = slower)",
            s.gridStep,
            0.06,
            0.25,
            (v) => setState(() => s = s.copyWith(gridStep: v)),
            (v) => v.toStringAsFixed(2),
          ),
          _slider(
            "Wall reflection",
            s.wallRefl,
            0.2,
            0.99,
            (v) => setState(() => s = s.copyWith(wallRefl: v)),
            (v) => v.toStringAsFixed(2),
          ),
          _slider(
            "Gain",
            s.gain,
            0.2,
            3.0,
            (v) => setState(() => s = s.copyWith(gain: v)),
            (v) => v.toStringAsFixed(2),
          ),
          _intSlider(
            "Reflections (0–4)",
            s.reflections,
            0,
            4,
            (v) => setState(() => s = s.copyWith(reflections: v)),
          ),
          const SizedBox(height: 10),
          SwitchListTile(
            contentPadding: EdgeInsets.zero,
            title: const Text("Robust heatmap scale"),
            subtitle: const Text(
              "Uses 98th percentile as 0 dB (reduces tiny-hotspot domination).",
            ),
            value: s.robustHeatmapScale,
            onChanged: (v) =>
                setState(() => s = s.copyWith(robustHeatmapScale: v)),
          ),

          const SizedBox(height: 12),
          TextField(
            controller: _optFreqCtrl,
            decoration: const InputDecoration(
              labelText: "Opt freqs (comma separated)",
              hintText: "e.g. 35,45,55,65",
            ),
          ),
          const SizedBox(height: 18),
          _sectionTitle(context, "Optimization"),
          _intSlider(
            "# speakers to optimize (1–6)",
            s.optSpeakers,
            1,
            6,
            (v) => setState(() => s = s.copyWith(optSpeakers: v)),
          ),
          _intSlider(
            "Iterations (more = slower)",
            s.optIters,
            100,
            4000,
            (v) => setState(() => s = s.copyWith(optIters: v)),
          ),

          const SizedBox(height: 10),
          DropdownButtonFormField<PlacementStyle>(
            initialValue: s.placementStyle,
            decoration: const InputDecoration(labelText: "Suggestion style"),
            items: const [
              DropdownMenuItem(
                value: PlacementStyle.realistic,
                child: Text("Realistic (recommended)"),
              ),
              DropdownMenuItem(
                value: PlacementStyle.free,
                child: Text("Free (experimental)"),
              ),
            ],
            onChanged: (v) {
              if (v == null) return;
              setState(() => s = s.copyWith(placementStyle: v));
            },
          ),
          const SizedBox(height: 18),
          _sectionTitle(context, "Seat"),
          DropdownButtonFormField<SeatType>(
            initialValue: s.seatType,
            decoration: const InputDecoration(labelText: "Seat type"),
            items: const [
              DropdownMenuItem(value: SeatType.couch, child: Text("Couch")),
              DropdownMenuItem(value: SeatType.chair, child: Text("Chair")),
            ],
            onChanged: (v) =>
                setState(() => s = s.copyWith(seatType: v ?? SeatType.couch)),
          ),

          const SizedBox(height: 18),
          _sectionTitle(context, "Room shape"),
          SwitchListTile(
            contentPadding: EdgeInsets.zero,
            title: const Text("Occlusion (heuristic)"),
            subtitle: const Text(
              "Attenuates speaker paths that leave the room outline. Helps non‑convex rooms (e.g. U‑shapes).",
            ),
            value: s.occlusionEnabled,
            onChanged: (v) =>
                setState(() => s = s.copyWith(occlusionEnabled: v)),
          ),
          SwitchListTile(
            contentPadding: EdgeInsets.zero,
            title: const Text("Soft corners (diffraction-ish)"),
            subtitle: const Text(
              "Reduces the occlusion penalty when the path passes near a corner (smoother shadows).",
            ),
            value: s.softDiffraction,
            onChanged: (v) =>
                setState(() => s = s.copyWith(softDiffraction: v)),
          ),
        ],
      ),
    );
  }

  Widget _sectionTitle(BuildContext context, String t) {
    return Padding(
      padding: const EdgeInsets.only(bottom: 8),
      child: Text(t, style: Theme.of(context).textTheme.titleMedium),
    );
  }

  Widget _slider(
    String label,
    double value,
    double min,
    double max,
    ValueChanged<double> onChanged,
    String Function(double) fmt,
  ) {
    return Padding(
      padding: const EdgeInsets.symmetric(vertical: 8),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text("$label: ${fmt(value)}"),
          Slider(value: value, min: min, max: max, onChanged: onChanged),
        ],
      ),
    );
  }

  Widget _intSlider(
    String label,
    int value,
    int min,
    int max,
    ValueChanged<int> onChanged,
  ) {
    return Padding(
      padding: const EdgeInsets.symmetric(vertical: 8),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text("$label: $value"),
          Slider(
            value: value.toDouble(),
            min: min.toDouble(),
            max: max.toDouble(),
            divisions: max - min,
            onChanged: (v) => onChanged(v.round()),
          ),
        ],
      ),
    );
  }
}

// ------------------------------------------------------------
// Helper struct
// ------------------------------------------------------------
class _ImageSource {
  final Offset pos;
  final int order;
  final double amp;
  const _ImageSource(this.pos, this.order, this.amp);
}
