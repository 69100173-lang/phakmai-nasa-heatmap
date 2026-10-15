// NASA/Google Earth Engine: Crop Water Stress Risk (prototype)
// วิธีใช้:
// 1) เปิด https://code.earthengine.google.com/ แล้วสร้าง Script ใหม่
// 2) วาด/นำเข้า polygon ขอบเขตพื้นที่ที่ได้รับการยืนยันจากหน่วยงาน
// 3) เปลี่ยนบรรทัด var area = geometry; ให้ตรงกับชื่อ geometry ของคุณ
// 4) Run แล้วตรวจชั้นข้อมูล ก่อน Export GeoTIFF ไป Google Drive
// ห้ามใช้กรอบพิกัดประมาณการแทนขอบเขตหมู่บ้านจริง
//
// ข้อจำกัด: เป็นดัชนีคัดกรองความเสี่ยงเบื้องต้น ไม่ใช่การวินิจฉัยพืชขาดน้ำ
// NDVI และ LST มีความละเอียดต่างกัน; ผลลัพธ์จะถูก resample ที่ 250 m
// ต้องตรวจเทียบกับข้อมูลภาคสนามก่อนใช้ตัดสินใจให้น้ำ

var area = geometry; // วาดหรือ import ขอบเขตที่ยืนยันแล้วใน GEE
var end = ee.Date(Date.now());
var start = end.advance(-32, 'day');
var baselineStart = end.advance(-5, 'year');

Map.centerObject(area, 13);
Map.addLayer(area, {color: 'white'}, 'ขอบเขตที่นำเข้า (ตรวจสอบแล้ว)');

function maskModisNdvi(img) {
  var qa = img.select('SummaryQA');
  return img.updateMask(qa.lte(1))
    .select('NDVI').multiply(0.0001).rename('ndvi')
    .copyProperties(img, ['system:time_start']);
}

var ndviCollection = ee.ImageCollection('MODIS/061/MOD13Q1')
  .filterBounds(area)
  .filterDate(baselineStart, end)
  .map(maskModisNdvi);

var ndviNow = ndviCollection.filterDate(start, end).median().select('ndvi');
var ndviBase = ndviCollection.filterDate(baselineStart, end.advance(-32, 'day'))
  .select('ndvi').mean();
var ndviDrop = ndviBase.subtract(ndviNow).divide(0.35).clamp(0, 1)
  .rename('ndvi_risk');

function maskLst(img) {
  var qa = img.select('QC_Day');
  // MOD11A2 QC_Day bits 0-1: mandatory QA; 0 means good quality
  var good = qa.bitwiseAnd(3).eq(0);
  return img.updateMask(good).select('LST_Day_1km')
    .multiply(0.02).subtract(273.15).rename('lst_c')
    .copyProperties(img, ['system:time_start']);
}

var lstCollection = ee.ImageCollection('MODIS/061/MOD11A2')
  .filterBounds(area).filterDate(baselineStart, end).map(maskLst);
var lstNow = lstCollection.filterDate(start, end).mean().select('lst_c');
var lstBase = lstCollection.filterDate(baselineStart, end.advance(-32, 'day'))
  .select('lst_c').mean();
var hotAnomaly = lstNow.subtract(lstBase).subtract(1).divide(8).clamp(0, 1)
  .rename('heat_risk');

// CHIRPS precipitation is a public rainfall estimate, not a NASA sensor.
// Compare the latest 30-day rainfall with the same calendar window in prior years.
var rainNow = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterBounds(area).filterDate(end.advance(-30, 'day'), end)
  .sum().rename('rain_mm');
var rainHistory = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterBounds(area)
  .filter(ee.Filter.calendarRange(end.get('month'), end.get('month'), 'month'))
  .filterDate(baselineStart, end.advance(-1, 'year'))
  .sum();
var rainBase = rainHistory.mean().rename('rain_baseline_mm');
var rainDeficit = rainBase.subtract(rainNow).divide(rainBase.max(1))
  .clamp(0, 1).rename('rain_risk');

// Initial prototype weights. Tune against field observations before operational use.
var risk = ndviDrop.multiply(0.35)
  .add(hotAnomaly.multiply(0.35))
  .add(rainDeficit.multiply(0.30))
  .multiply(100).clamp(0, 100).rename('water_stress_risk');

// MODIS LST is coarser than NDVI. Export at 250 m for a consistent raster grid;
// this does not create 250 m thermal detail.
var risk250 = risk
  .resample('bilinear')
  .clip(area)
  .toFloat();

var palette = ['218c4a', 'e5d83b', 'ef982f', 'd93636'];
Map.addLayer(ndviNow.clip(area), {min: 0, max: 0.9, palette: ['8c510a','f6e8c3','5ab4ac','01665e']}, 'NDVI ล่าสุด');
Map.addLayer(lstNow.clip(area), {min: 20, max: 45, palette: ['2c7bb6','abd9e9','ffffbf','fdae61','d7191c']}, 'LST กลางวัน °C');
Map.addLayer(risk250, {min: 0, max: 100, palette: palette}, 'ความเสี่ยงขาดน้ำ 0–100', true, 0.75);

print('ช่วงข้อมูลล่าสุด', start, end);
print('จำนวนภาพ NDVI ในช่วงล่าสุด', ndviCollection.filterDate(start, end).size());
print('จำนวนภาพ LST ในช่วงล่าสุด', lstCollection.filterDate(start, end).size());
print('คะแนนความเสี่ยง (สถิติในพื้นที่)', risk250.reduceRegion({
  reducer: ee.Reducer.minMax().combine(ee.Reducer.mean(), '', true),
  geometry: area, scale: 1000, maxPixels: 1e8
}));

Export.image.toDrive({
  image: risk250,
  description: 'Phakmai_Water_Stress_Risk',
  folder: 'EarthEngineExports',
  fileNamePrefix: 'phakmai-water-stress-risk',
  region: area,
  scale: 250,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});
