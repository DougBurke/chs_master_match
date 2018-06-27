/*
 * Code for the master-hull page.
 */

"use strict";

var settings;

var state;

var pageNum;
var imageScale = 'log10';

const psfColor = 'white';

const defBinACIS = 8;
const defBinHRC = 32;

// Change to the given page and scaling.
//
function changePage() {
  let pageStr = pageNum.toString();
  if (pageNum < 10) { pageStr = "00" + pageStr; }
  else if (pageNum < 100) { pageStr = "0" + pageStr; }

  const img = '/img/' + settings.ensemble + '/hull.' +
      settings.ensemble + '.' +
      settings.mid + '.p' + pageStr +
      '.v' + settings.revstr + '.' +
      imageScale + '.png';

  document.getElementById('pageview').src = img;
}

function setPage(newPage) {
  if ((newPage < 1) || (newPage > settings.npages)) {
    alert("Internal error: newPage = " + newPage.toString());
    return;
  }

  pageNum = newPage;
  changePage();
}

function setScaling(newScale) {
  imageScale = newScale;
  changePage();
}


// See https://stackoverflow.com/a/10284006
//
function zip() {
  var args = [].slice.call(arguments);
  var shortest = args.length==0 ? [] : args.reduce((a,b) => {
    return a.length < b.length ? a : b
  });

  return shortest.map(function(_,i){
    return args.map((array) => {return array[i];});
  });
}

// Store the current coordinates of the master hulls - i.e. the
// user selection before saving. This is updatePage by
// initialize and modified by handleRegionChange().
//
var masterhulls_raw = [];
var masterhulls_cnvx = [];

// Need mapping from which hull for this master and the polygon id
// (for the QA case where there are multiple master hulls), and
// the polygon id to this index (so that we can change things between
// different JS9 windows). This is actually overly complex, and just guards
// against issues where the user may have added regions.
//
var masterhulls_ids = {
    // store the polygon id for each stack/cpt key
    master: {}, convex: {},
    // store the index for each polygon id
    master_index: {}, convex_index: {}
};

// pan to this position; opts is the argument to pass
// to JS9 commands to determine the window to use.
// wcs contains ra and dec fields.
//
function goToRaDec(wcs, opts, debug=true) {
  const pix = JS9.WCSToPix(wcs, opts);
  if (debug) {
      console.log("-> " + wcs.ra + " " + wcs.dec + " : " + pix.x + " " + pix.y);
  }
  JS9.SetPan(pix, opts);
}

// This is a wrapper around JS9.DisplaySection which restores the
// current pan position after a change. Hopefully (as of JS9 v2.1)
// this should not be needed, but I have not validated this assertion.
//
function changeJS9Display(stack, cpt, section, opts) {

  // What is the location of the center?
  const pix = JS9.GetPan(opts);

  // Convert to WCS in case the binning is changed
  const wcs = JS9.PixToWCS(pix, opts);

  const key = toKey(stack, cpt);
  const blur = blurVal[key];

  const args = Object.assign({}, section);
  args.ondisplaysection = (im) => {
      const newpix = JS9.WCSToPix(wcs, opts);
      JS9.SetPan(newpix, opts);
      if (typeof blur !== "undefined" ) {
        JS9.GaussBlurData(blur, opts);
      }
  };
  JS9.DisplaySection(args, opts);

}

// Finalize the widgets in the window and add the regions.
//
function finalizeJS9Display(stack, cptnum, winid) {

  return function(img) {
    setupJS9(img, stack, cptnum, winid);
    addRegionsToJS9(img, stack, cptnum, winid);
  };
}

//   - psf        for PSF ellipses
//   - stack      for stack-level hulls
//   - original   the input 'master-level' hulls
//   - convex     the convex hull calculated from the current hull
//   - regions    the user-editable hull
//                i.e. this is the default regions layer
//
// Only the latter is user-editable or selectable (at least
// directly). For now we use the "standard" regions layer
// for this
//
const psfLayer = 'psf';
const stackLayer = 'stack';
const convexLayer = 'convex';
const originalLayer = 'original';
const masterLayer = 'regions';

// Need a dictionary key for stack,cpt, so convert to
// a string.
//
function toKey(stack, cpt) {
  return  stack + "." + cpt.toString();
}

// Change the blurring for the data displayed in the JS9 instance.
//
// We record the value of the blurring so that it can be restored
// when the bin/band/... are changed, but this relies on the blur value
// *only* being changed by our button and not the JS9 menu item.
//
let blurVal = {};
function blurData(stack, cpt, newval, opts) {
  const key = toKey(stack, cpt);
  blurVal[key] = newval;
  JS9.GaussBlurData(newval, opts);
}

// Change the binning for the data displayed in the JS9 instance.
//
function binData(stack, cpt, newval, opts) {
  changeJS9Display(stack, cpt, {bin: newval}, opts);
}

// Customize the JS9 display window
//
// winid is the base HTML id of the "light window" containing
// the JS9 display.
//
// Is also sets up the layers.
//
function setupJS9(img, stack, cptnum, winid) {

  // Hard code logic for the parent element (could find it by
  // hunting up the tree but can not be bothered at the moment)
  //
  const container = document.getElementById('d' + winid);

  // Add handlers for the user buttons:
  //   - blur
  //   - rebin
  //   - save region
  //   - reload region
  //   - zoom
  //
  const opts = {display: img};
  let btns = container.getElementsByClassName("sigma");
  for (let i = 0; i < btns.length; i++) {
    btns[i].addEventListener("change", (e) => {
      blurData(stack, cptnum, e.target.value, opts);
    });
  }

  btns = container.getElementsByClassName("binsize");
  for (let i = 0; i < btns.length; i++) {
    btns[i].addEventListener("change", (e) => {
      binData(stack, cptnum, e.target.value, opts);
    });
  }

  // Save the current version of the hull.
  //
  const saveButton = document.getElementById(winid + 'SaveMasters');
  if (saveButton !== null) {
    saveButton
      .addEventListener("click",
			(e) => { saveMasters(); });
  }

  // Go back to the last-saved version of the master hull.
  //
  // I do not want to delete existing regions since I want the
  // displays to automatically change - and rely on the change
  // handler to do this - *but* this only works for the "proposed"
  // shape, and not the "original" or "convex hull" shapes.
  //
  // However, the "original" layer should not be relevant here,
  // since it should not be changed by this call (it's the proposed
  // shape, which hasn't been changed by the user).
  //
  // The convex hull is based on the hull layer, so it's not
  // clear how this is picked up/handled.
  //

  const reloadButton = document.getElementById(winid + 'ReloadMasters');
  if (reloadButton !== null) {
    reloadButton
      .addEventListener("click",
			(e) => {
			  setupMasters();
			  addMasterHullToJS9(opts, stack, cptnum, true);
			});
  }

  const bandSelect = document.getElementById(winid + 'BandChoice');
  if (bandSelect !== null) {
      bandSelect
	  .addEventListener("change",
			    (e) => { changeFilter(stack, cptnum,
						  e.target.value, opts);
			    });
  }

  const psfColSelect = document.getElementById(winid + 'PSFColor');
  if (psfColSelect !== null) {
      psfColSelect
	  .addEventListener("change", (e) => {
		  colorizePSFs(stack, winid, e.target.value);
	      });
  }

  container.getElementsByClassName("zoomin")[0]
     .addEventListener("click", (e) => { JS9.SetZoom("in", opts); });
  container.getElementsByClassName("zoomout")[0]
     .addEventListener("click", (e) => { JS9.SetZoom("out", opts); });

  // toggle the panner; it looks like need to give winid to display,
  // not img.
  //
  document.getElementById(winid + 'ShowPanner')
      .addEventListener("click",
			(e) => { JS9.DisplayPlugin('panner',
                                                   {display: winid}); });

  // Set up the layers. We have to base the options on
  // JS9.Regions.opts to get "sensible" behavior.
  //
  let layerOpts = Object.assign({}, JS9.Regions.opts);
  for (const name of ['movable', 'rotatable', 'resizable', 'evented']) {
    layerOpts[name] = false;
  }

  // Make sure we remove the onchange handler from these extra layers
  layerOpts.onchange = null;

  // Eric hopes this will address the half-pixel offset on rebin
  // but I need to investigate it further to find out what is
  // going in.
  //
  // layerOpts.dowcsstr = true;

  for (const name of [psfLayer, stackLayer, originalLayer, convexLayer]) {
    JS9.NewShapeLayer(name, layerOpts, opts);
  }
}

// I'm adding in regions one at a time; the JS9 API lets you specify
// multiple regions at once, so perhaps I would be better off
// following that approach - i.e. create all the shapes for a
// given layer and then add or change them. However, I end up
// storing the id value returned by AddShapes, to use later,
// so this is not possible, at least for those that I want to
// change.
//
function add_hull_to_js9(hull, opts, win, layer='regions') {

  // Need to convert to image coordinates
  const pts = [];
  for (let i = 0; i < hull.ra.length; i++) {
    pts.push(JS9.WCSToPix(hull.ra[i], hull.dec[i], win));
  }

  // Need to set the label
  const shape = {shape: 'polygon', pts: pts};
  if (typeof(hull.label) !== "undefined") {
    shape.text = hull.label;
  }
  return JS9.AddShapes(layer, shape, opts, win);

}


function change_hull_in_js9(hull, layer, idval, win) {

  // Need to convert to image coordinates
  const pts = [];
  for (let i = 0; i < hull.ra.length; i++) {
    pts.push(JS9.WCSToPix(hull.ra[i], hull.dec[i], win));
  }
  JS9.ChangeShapes(layer, idval, {pts: pts}, win);
}


function addPSFRegions(stack, win) {
  const psfs = settings.regionstore.stackpsfs[stack];
  if (typeof psfs === "undefined") {
    return;
  }

  const psfOpts = {color: psfColor,
		   changeable: false,
		   tags: 'psf'};

  // Let JS9 deal with the WCS conversion rather than doing it
  // here. My guess is that the region-parsing that JS9 does is
  // going to be as efficient as me doing it manually (or at least
  // not significantly different in terms of the runtime).
  //
  for (let psf of psfs) {
    const shape = 'fk5; ellipse(' + psf.ra.toString() + ',' + 
	psf.dec.toString() + ',' +
        psf.r0.toString() + '",' +
        psf.r1.toString() + '",' +
        psf.angle.toString() + ')';

    JS9.AddShapes(psfLayer, shape, psfOpts, win);
  }
}

// Note that newcol can be set to 'hide', which means toggle off.
//
function colorizePSFs(stack, winid, newcol) {
  if (newcol === 'hide') {
    JS9.ShowShapeLayer(psfLayer, false, {display: winid});
  } else {
    // could track whether a show call is needed, but be lazy
    JS9.ShowShapeLayer(psfLayer, true, {display: winid});
    JS9.ChangeShapes(psfLayer, "all", {color: newcol}, {display: winid});
  }
}

// Change the filter to the given band.
//
// The stack and cpt are sent in to allow the current
// position to be retained, if desired.
//
function changeFilter(stack, cpt, band, opts) {
  const enfilter = band_to_filter(band);
  changeJS9Display(stack, cpt, {filter: enfilter}, opts);
}

function makeConvexOpts(name) {
  return {color: 'lime',
	  strokeDashArray: [5, 1],
	  changeable: false,
	  evented: false,
	  tags: name};
}

function makeMasterOpts(name, isqa) {
  var color;
  if (isqa) { color = 'cyan'; } else { color = 'gold'; }

  return {movable: false,
	  color: color,
	  rotatable: false,
	  resizable: false,
	  tags: name};
}

// This function adds regions (when reload is false) or changes
// regions (when reload is true).
//
// TODO: Deleting stuff is also going to mess things up.
//
function addMasterHullToJS9(display, stack, cptnum, reload=false) {

console.log("In addMasterHullToJS9");

  const masterhull = settings.regionstore.masterhulls[settings.masterid];
  if (typeof masterhull === "undefined") {
    alert("No master hull found!");
    return;
  }

  if (masterhull.wcs.length === 0) {
    console.log("We have no master hull information for " +
                settings.ensemble + " " +
		settings.masterid.toString() + "!");
    return;
  }

  const is_qa = masterhull.status.startsWith('qa');

  if (reload) {
    redrawMasterHulls(display, stack, cptnum);
  } else {
      addMasterHulls(display, stack, cptnum, is_qa);
  }
}

function redrawMasterHulls(display, stack, cptnum) {

  const key = toKey(stack, cptnum);
  let idstore = masterhulls_ids.convex[key];
  for (var i = 0; i < masterhulls_cnvx.length; i++) {
    change_hull_in_js9(masterhulls_cnvx[i], convexLayer, idstore[i], display);
  }

  idstore = masterhulls_ids.master[key];
  for (var i = 0; i < masterhulls_raw.length; i++) {
    change_hull_in_js9(masterhulls_raw[i], masterLayer, idstore[i], display);
  }
}

function addMasterHulls(display, stack, cptnum, is_qa) {

  const key = toKey(stack, cptnum);
  if ((key in masterhulls_ids.convex) || (key in masterhulls_ids.master)) {
    alert("Internal error: key=" + key + " is already in masterhulls_ids");
  }
  
  // Note: the tagName gets over-written later, so should change this.
  const tagName = 'master';
  const hullOpts = makeMasterOpts(tagName, is_qa);

  /* If this is a review page then don't let the user change the hull */
  if (state.ensemble_status !== "todo") {
    hullOpts.changeable = false;
  }

  // Note: the name will be changed later.
  const convexName = 'master';
  const convexOpts = makeConvexOpts(convexName);

  masterhulls_ids.convex[key] = [];
  masterhulls_ids.convex_index[key] = [];
  for (var i = 0; i < masterhulls_cnvx.length; i++) {
    const hull = masterhulls_cnvx[i];

    const shapeName = convexName + i.toString();
    convexOpts.tags = [shapeName, convexName];

    const polyid = add_hull_to_js9(hull, convexOpts, display, convexLayer);
    masterhulls_ids.convex[key][i] = polyid;
    masterhulls_ids.convex_index[key][polyid] = i;
  }

  masterhulls_ids.master[key] = [];
  masterhulls_ids.master_index[key] = [];
  for (var i = 0; i < masterhulls_raw.length; i++) {
    const hull = masterhulls_raw[i];

    const shapeName = tagName + i.toString();
    hullOpts.tags = [shapeName, tagName];

    const polyid = add_hull_to_js9(hull, hullOpts, display, masterLayer);
    masterhulls_ids.master[key][i] = polyid;
    masterhulls_ids.master_index[key][polyid] = i;
  }

}

// Add the stack-level and master hull(s) to the JS9 window.
//
// The PSF regions are drawn first (if available).
//
// Stack level hulls are drawn first (maybe in a different layer?)
//
function addRegionsToJS9(img, stack, cptnum, regions) {

  let stackhulls = settings.regionstore.stackhulls[stack];
  if (typeof stackhulls === "undefined") {
    alert("no stack-level hulls found for " + stack);
    return;
  }

  const display = {display: img};

  addPSFRegions(stack, display);

  let linestyle = [1];
  const hullOpts = {color: 'orange',
		    strokeDashArray: linestyle,
		    changeable: false,
		    tags: 'stack'};

  let ra0 = undefined, dec0 = undefined;
  for (const shull of stackhulls) {

    // linestyle: solid for mancode is 0, otherwise
    // dotted.
    //
    // wanted to add strokeWidth to increase the width, but can
    // not get it to work
    //
    // The color depends on whether this is the "selected" component
    // or not.
    //
    if (shull.component == cptnum) {
      hullOpts.color = 'orange';
      ra0 = shull.ra0;
      dec0 = shull.dec0;
    } else {
      hullOpts.color = '#cc3333'; // a red-ish color
    }

    if (shull.mancode === 0) {
      delete hullOpts.strokeDashArray;
    } else {
      hullOpts.strokeDashArray = linestyle;
    }

    add_hull_to_js9(shull, hullOpts, display, stackLayer);
  }

  // The "original" shape can not change, and so we draw it here
  // rather than addMasterHullToJS9, which deals with those things
  // that can change (either because the user has tweaked the
  // settings or reloaded them from the last-saved values).
  //
  const masterhull = settings.regionstore.masterhulls[settings.masterid];

  const tagName = "original";
  const origOpts = makeMasterOpts(tagName);
  origOpts.color = 'white';
  origOpts.strokeDashArray = [3, 2];

  // NOTE use of wcs_orig so we always see the proposed hull as
  //      the original one
  for (var i = 0; i < masterhull.wcs_orig.length; i++) {
      const hull = masterhull.wcs_orig[i];

      // Setting the name here is suboptimal
      const shapeName = tagName + i.toString();
      origOpts.tags = [shapeName, tagName];

      add_hull_to_js9(hull, origOpts, display, originalLayer);
  }

  addMasterHullToJS9(display, stack, cptnum);

  // Move to the center of the stack-level hull, if defined,
  // rather than the center of the master hull.
  //
  if (typeof ra0 !== "undefined") {
    goToRaDec({ra: ra0, dec: dec0}, display);
  }
}

function blur_html(id) {
  const name = id + "sigma";
  let html = "<div class='blur'><span>Blur</span>";
  const sigmas = [0, 1, 2, 3, 4];
  for (let i = 0; i < sigmas.length; i++) {
    /* need to have label after the input for the CSS */
    const sigma = sigmas[i].toString();
    var l = id + 'sigma' + sigma;
    html += "<input class='sigma' name='" + name + "' type='radio' ";
    html += "id='" + l + "' value='" + sigma + "'"
    if (i === 0) { html += " checked"; }
    html += ">";
    html += "<label for='" + l + "'>" + sigma + "</label>";
  }
  html += "</div>";
  return html;
}

function bin_html(stack, id) {
  const all_bins = [1, 2, 4, 8, 16, 32, 64, 128];
  let def_binsize, def_start;
  if (stack.startsWith('hrc')) {
    def_binsize = defBinHRC;
    def_start = 1;
  } else {
    def_binsize = defBinACIS;
    def_start = 0;
  }
  const def_bins = all_bins.slice(def_start, def_start + 7);

  const name = id + "binsize";
  let html = "<div class='rebin'><span>Bin</span>";
  for (let i = 0; i < def_bins.length; i++) {
    const binsize = def_bins[i].toString();
    const l = id + 'binsize' + binsize;
    html += "<input class='binsize' name='" + name + "' type='radio' ";
    html += "id='" + l + "' value='" + binsize + "'";
    if (def_bins[i] === def_binsize) {
      html += " checked";
    }

    html += ">";
    html += "<label for='" + l + "'>" + binsize + "</label>";
  }
  html += "</div>";
  return html;
}

function zoom_html() {
  let html = "<span class='zoom'>";
  html += "<button class='zoomin'>+</button>";
  html += "<button class='zoomout'>-</button>";
  html += "</span>";
  return html;
}

function panner_html(id) {
  let html = "<span class='panner'>"; // "Toggle: ";
  html += "<button id='" + id + "ShowPanner'>Panner</button>";
  html += "</span>";
  return html;
}

function psf_html(stack, id) {
  const psfs = settings.regionstore.stackpsfs[stack];
  if (typeof psfs === "undefined") { return ""; }

  let html = "<span class='colorize'>PSF: ";
  html += "<select id='" + id + "PSFColor'>";
  for (const newopt of ['hide', 'white', 'yellow', 'black', 'red',
  			    'orange', 'cyan', 'blue', 'brown']) {
  	  html += "<option value='" + newopt + "'";
  	  if (newopt === psfColor) {
  	      html += " selected";
  	  }
  	  html += ">" + newopt + "</option>";
  }

  html += "</select>";
  html += "</span>";

  return html;
}

function band_html(band, id) {
  if (band === "w") { return ""; }

  let html = "<span class='bandchoice'>Band: ";
  html += "<select id='" + id + "BandChoice'>";

  for (const bandval of ['b', 'u', 's', 'm', 'h']) {
  	  html += "<option value='" + bandval + "'";
  	  if (bandval === band) {
  	      html += " selected";
  	  }
  	  html += ">" + bandval + "</option>";
  }

  html += "</select>";
  html += "</span>";
  
  return html;
}

function save_html(id) {
  let html = "<span class='save'>";
  html += "<button id='" + id;
  html += "SaveMasters'>Save</button>";
  html += "</span>";
  return html;
}

function reload_html(id) {
  let html = "<span class='reload'>";
  html += "<button id='" + id;
  html += "ReloadMasters'>Reload</button>";
  html += "</span>";
  return html;
}

/*
 * HTML code for the JS9 display.
 *
 * Note that the event listeners are added in the
 * setupJS9 callback.
 */
function js9_display_html(stack, stacknum, cptnum, band, id) {
  let html = "<div class='stackid' id='" + id + "stackid'>" +
    "<span style='float: left;'>Stack: " + 
      stacknum.toString() + ": " + stack +
      " cpt: " + cptnum.toString() + 
      " band: " + band +
      "</span></div>";

  html += "<div class='useropts'>";

  html += "<div class='kitchensink'>";

  // Have the "always-have" options first, and the optional
  // ones last.
  html += save_html(id);
  html += reload_html(id);

  html += band_html(band, id);
  html += psf_html(stack, id);

  html += "</div>"; // class=kitchensink

  html += "<div class='buttonopts'>";
  html += blur_html(id);
  html += bin_html(stack, id);
  html += zoom_html();
  html += panner_html(id);
  html += "</div>";

  html += "</div>"; // class=useropts

  html += "<div class='JS9Menubar' id='" + id + "Menubar'></div>";
  html += "<div class='JS9' id='" + id + "'></div>'";
  html += "<div class='JS9Colorbar' id='" + id + "'></div>'";

  return html;
}

// Store the mapping between the window id and the stack + component.
//
var js9_mapping = {};

// "unique" id for new JS9 windows
//
//
var idctr = 1;
function js9_id(stack, cptnum) {

  var retval = "js9win" + idctr.toString();
  if (retval in js9_mapping) {
      alert("Internal error: " + retval + " is not unique!");
  }
  js9_mapping[retval] = {stack: stack, cpt: cptnum};

  idctr += 1;
  return retval;
}

// Create the filter expression for the band.
//
// Is the minimum energy for the ultra-soft band correct?
// I don't think we have any u-band CHS so it is somewhat academic.
//
function band_to_filter(band) {
    if (band == 'b')      { return "energy >= 500 && energy < 7000"; }
    else if (band == 'u') { return "energy >= 200 && energy < 500"; }
    else if (band == 's') { return "energy >= 500 && energy < 1200"; }
    else if (band == 'm') { return "energy >= 1200 && energy < 2000"; }
    else if (band == 'h') { return "energy >= 2000 && energy < 7000"; }
    else if (band == 'w') { return ""; }
    else {
	console.log("Unexpected band: " + band);
	return "";
    }
}

// Applies an energy filter for ACIS data to try and match the hull(s).
//

function showInJS9(val) {
  if (val.trim() === '') { return; }

  const toks = val.split(',');
  if (toks.length != 3) {
    alert("Internal error: val=[" + val + "]");
    return;
  }

  const stack = toks[0];
  const stacknum = toks[1];
  const cptnum = toks[2];

  const winid = js9_id(stack, cptnum);
  const opts = {onload: finalizeJS9Display(stack, cptnum, winid)};

  const band = settings.enbands_cpt[stack + "." + cptnum];

  if (stack.startsWith('acis')) {
    opts.bin = defBinACIS;
    opts.xdim = 8192;
    opts.ydim = 8192;
    // does setting this actually do anything?
    opts.filter = band_to_filter(band);
  } else {
    opts.bin = defBinHRC;
    opts.xdim = 16384;
    opts.ydim = 16384;
  }
  opts.id = winid;
  return JS9.LoadWindow('/evt3/' + stack, opts, 'light',
                        js9_display_html(stack, stacknum, cptnum, band,
					 opts.id));
}

/*** callback code for handling region changes ***/

/*
 * Let's see how well it works with the onchange handler.
 * If not, could we chain on the mouseup handler?
 */

// The code adds tags for master and stack which could be used
// to distinguish regions, although stack polygons should not
// be changeable.
//
// What about regions that the user has added directly?
//
function handleRegionChange(img, action) {
  console.log(">> mode=" + action.mode + " tags=" + action.tags);

  // For now, we only care about master tags;
  // this may need to be tweaked if we use the tag name as an id
  //
  if (!action.tags.includes('master')) { return; }

  if (action.mode === "update") {
    broadcastMasterUpdate(img, action);
  } else if (action.mode === "remove") {
    broadcastMasterDelete(img, action);
  }
}

// Tell the other windows about the new polygon.
//
function broadcastMasterUpdate(img, action) {

  // This is the window that the user has changed, so we don't want
  // to change this window.
  const js9win = img.display.id;
  const baseWin = {display: js9win};

  // What is the "index" for this polygon?
  const baseHullinfo = js9_mapping[js9win];
  const baseKey = toKey(baseHullinfo.stack, baseHullinfo.cpt);

  // These two should be the same
  const baseConvexIndex = masterhulls_ids.convex_index[baseKey][action.id];
  const baseMasterIndex = masterhulls_ids.master_index[baseKey][action.id];

  // console.log("indexes into array: " + baseConvexIndex + " " + baseMasterIndex);

  // Convert the polygon into a convex hull (if necessary).
  // I want to do this conversion in a "tangent-plane" projection
  // where I don't have to deal with ra=0/360 wrap around.
  //
  let chull_sky = window.convexHull(action.pts);
  let chull_eqpos = [];
  for (let sky of chull_sky) {
    chull_eqpos.push(JS9.PixToWCS(sky, baseWin));
  }

  // Since using LightWindows, can look for div.dhtmlwindow
  // containers, and the knowledge that the id for this is
  // 'd' + id-of-js9-div
  //
  // Although this doesn't work as soon as the user creates a panner
  // or other related window. These extra divs appear to have
  // ids like js9win2_JS9Panner_lightDiv
  // as opposed to djs9win2, which we want. So for now just
  // reject anything with an underscore.
  //
  // The convex hull is
  //   a) drawn in all windows, including the one being edited
  //   b) drawn first, so that it hopefully appears below the
  //      master hull, so that the master hull can still be
  //      selected/edited.
  //

  const divs = document.getElementsByClassName('dhtmlwindow');
  for (let div of divs) {
    if (div.id.includes('_')) {
      continue;
    }
    const owin = div.id.substring(1);

    const imname = {display: owin};
    const hdl = JS9.GetImage(imname);

    // What stack/cpt is this?
    const hullinfo = js9_mapping[owin];
    const key = toKey(hullinfo.stack, hullinfo.cpt);

    // What region ids do we want?
    const convexId = masterhulls_ids.convex[key][baseConvexIndex];
    const masterId = masterhulls_ids.master[key][baseMasterIndex];

    // stop propogating the onchange signal
    hdl.params.xeqonchange = !hdl.params.xeqonchange;

    const hullpts = [];
    for (const wcs of chull_eqpos) {
      hullpts.push(JS9.WCSToPix(wcs, imname));
    }

    JS9.ChangeShapes(convexLayer, convexId, {pts: hullpts}, imname);

    // Only adjust the master polygon if this is not the
    // window the user is changing.
    if (owin !== js9win) {

      // Need to convert to image coordinates
      const pts = [];
      for (const wcs of action.wcspts) {
        pts.push(JS9.WCSToPix(wcs, imname));
      }

      JS9.ChangeShapes(masterLayer, masterId, {pts: pts}, imname);
    }

    hdl.params.xeqonchange = !hdl.params.xeqonchange;

    // Ensure the regions layer is where the user action happens
    JS9.ActiveShapeLayer(masterLayer, imname);
  }

  // Store the new polygon (in Equatorial coordinates) in the
  // masterhulls_raw and masterhulls_cnvx global variables.
  //
  masterhulls_raw[baseMasterIndex] = {ra: [], dec: []};
  for (const wcs of action.wcspts) {
    masterhulls_raw[baseMasterIndex].ra.push(wcs.ra);
    masterhulls_raw[baseMasterIndex].dec.push(wcs.dec);
  }

  masterhulls_cnvx[baseConvexIndex] = {ra: [], dec: []};
  for (const wcs of chull_eqpos) {
    masterhulls_cnvx[baseConvexIndex].ra.push(wcs.ra);
    masterhulls_cnvx[baseConvexIndex].dec.push(wcs.dec);
  }

}

// TODO: this copies from broadcastMasterUpdate, so need
//       to refactor
//
// Deletes the master shape from both the region and
// convex layers. It does not delete the original version.
//
function broadcastMasterDelete(img, action) {

  const js9win = img.display.id;

  // Since using LightWindows, can look for div.dhtmlwindow
  // containers, and the knowledge that the id for this is
  // 'd' + id-of-js9-div
  //
  // Although this doesn't work as soon as the user creates a panner
  // or other related window. These extra divs appear to have
  // ids like js9win2_JS9Panner_lightDiv
  // as ooposed to djs9win2, which we want. So for now just
  // reject anything with an underscore.
  //
  const divs = document.getElementsByClassName('dhtmlwindow');
  for (let div of divs) {
    if (div.id.includes('_')) {
      continue;
    }
    const owin = div.id.substring(1);
    if (owin === js9win) {
      continue;
    }

    const imname = {display: owin};
    let hdl = JS9.GetImage(imname);

    // stop propogating the onchange signal
    hdl.params.xeqonchange = !hdl.params.xeqonchange;

    // JS9.RemoveRegions('master', imname);
    for (let layer of [convexLayer, masterLayer]) {
      JS9.RemoveShapes(layer, 'master', imname);
    }

    hdl.params.xeqonchange = !hdl.params.xeqonchange;
  }

  // Delete the convex-hull version in the original window.
  JS9.RemoveShapes(convexLayer, 'master', {display: js9win});

}

function saveMasters() {
  let httpRequest = new XMLHttpRequest();
  if (!httpRequest) {
      alert("Unable to create a XMLHttpRequest!");
      return;
  }

  // for now just store the hull itself, and not the convex hull
  // version.
  //
  const store = {ensemble: settings.ensemble,
		 revision: settings.revstr,
		 masterid: settings.masterid,
		 polygons: masterhulls_raw};

  // Add the spinner to the whole page
  //
  let body = document.getElementsByTagName("body")[0];
  let spinner = new Spinner(spinopts);
  spinner.spin(body);

  httpRequest.addEventListener("load", function() {
  });
  httpRequest.addEventListener("error", function() {
      alert("Unable to save data!");
  });
  httpRequest.addEventListener("abort", function() {
      alert("Unable to save data!");
  });
  httpRequest.addEventListener("loadend", function() {
      spinner.stop();
  });

  httpRequest.open('POST', '/save/masterpoly');
  httpRequest.setRequestHeader("Content-Type", "application/json;charset=UTF-8");
  httpRequest.send(JSON.stringify(store));
}

// Save user selections
//
function saveUser(store, field, newval) {
  let httpRequest = new XMLHttpRequest();
  if (!httpRequest) {
      alert("Unable to create a XMLHttpRequest!");
      return;
  }

  // Add the spinner to the whole page
  //
  let body = document.getElementsByTagName("body")[0];
  let spinner = new Spinner(spinopts);
  spinner.spin(body);

  // Update the local state once the save has been made. Note this happens
  // even during an error, I think.
  // Since both fields can be changed, change both of them.
  //
  httpRequest.addEventListener("load", function() {
    state.useraction.user = store.useraction;
    state.usernotes.user = store.usernotes;
  });
  httpRequest.addEventListener("error", function() {
      alert("Unable to save data!");
  });
  httpRequest.addEventListener("abort", function() {
      alert("Unable to save data!");
  });
  httpRequest.addEventListener("loadend", function() {
      spinner.stop();
  });

  httpRequest.open('POST', '/save/master');
  httpRequest.setRequestHeader("Content-Type", "application/json;charset=UTF-8");
  httpRequest.send(JSON.stringify(store));
  
}

function saveUserContent() {
  const newval = document.getElementById('usercontent').value;
  if (newval === state.usernotes.user) {
    return;
  }
  const store = {ensemble: settings.ensemble,
		 revision: settings.revstr,
		 masterid: settings.masterid,
		 useraction: "",
		 usernotes: newval};

  // You should only be able to save concent if this is the
  // latest version, so choice should always be set here.
  //
  const choice = document.getElementById('choice');
  if (choice !== null) {
      store.useraction = choice.value;
  }
  saveUser(store, "usernotes", newval);
}

// Note: this also saves any changes to the user content.
//       It could ignore that, but that seems surprising to the
//       user, and makes the "choice" button a "final arbiter".
//
function saveUserChoice(newval) {
  const notes = document.getElementById('usercontent').value;
  const store = {ensemble: settings.ensemble,
		 revision: settings.revstr,
		 masterid: settings.masterid,
		 useraction: newval,
		 usernotes: notes};
  saveUser(store);
}

/* play with websamp */

var sampConnection;

// UGLY
var sampBaseUrl = window.location.href.toString()
                     .replace(new RegExp("/[^/]*/[^/]*/[^/]*$"), "");

function showInSAMP(stack) {
  if (stack.trim() === '') { return; }
  const url = sampBaseUrl + "/evt3/" + stack;
  console.log("*** SAMP : " + url);

  sampConnection.runWithConnection((connection) => {
    var msg = new samp.Message("image.load.fits",
                               {"url": url, "name": stack});
    connection.notifyAll([msg]);
  });
}

function sampIsAvailable(flag) {
  document.getElementById("sampcontrols").hidden = !flag;
}

function finalize() {
  if (typeof sampConnection !== "undefined") {
    sampConnection.unregister();
    sampConnection = undefined;
  }
}

// Copy the "saved" master hulls into the "current" master hulls.
//
function setupMasters() {

  // Store away the current master hull settings for use by
  // addMasterHullToJS9 and handleRegionChange.
  //
  const masterhull = settings.regionstore.masterhulls[settings.masterid];
  masterhulls_raw = [];  // should not be needed
  for (const hull of masterhull.wcs) {
    masterhulls_raw.push(Object.assign({}, hull));
  }

  // I would like to assume that the hull is convex, so the coordinates
  // could be copied into masterhulls_cnvx, but I don't want to enforce
  // this until the user explicitly requests it. So, the masterhull
  // may not be convex.
  //
  // The problem is that there is no image to convert the RA/Dec
  // values to SKY, so I do not really want to convert to a convex
  // hull here. Although perhaps I should since it shouldn't really
  // matter too much (famous last words). The intention is that the
  // convex hull will get recalculated before use, but we want
  // something here just in case. Whether that is true is something
  // we shall find out.
  //
  masterhulls_cnvx = [];
  for (const hull of masterhull.wcs) {
    const inpts = [];
    for (let i = 0; i < hull.ra.length; i++) {
      inpts.push({x: hull.ra[i], y: hull.dec[i]});
    }
    const outpts = window.convexHull(inpts);
    const out = {ra: [], dec: []};
    for (let i = 0; i < outpts.length; i++) {
      out.ra.push(outpts[i].x);
      out.dec.push(outpts[i].y);
    }
    masterhulls_cnvx.push(out);
  }

}

var colorsShown = false;
function toggleColorMapping() {
  const el = document.getElementById('colormapping');
  var lbl;
  if (colorsShown) {
    el.style.display = "none";
    lbl = 'Show colors';
  } else {
    el.style.display = "block";
    lbl = 'Hide colors';
  }
  colorsShown = !colorsShown;
  document.getElementById('colormappingselect').innerHTML = lbl;
}

// Set up the page

function updatePage(json) {

  state = Object.assign({}, json);

  setupMasters();

  // Do we need a handler for region changes?
  if (state.ensemble_status === "todo") {
    if (JS9.Regions.opts.onchange !== null) {
      alert("Overwriting onchange handler!");
    }
    JS9.Regions.opts.onchange = handleRegionChange;
  }

  document.getElementById('imgscale').
    addEventListener("change", (e) => { setScaling(e.target.value); });

  document.getElementById('showinsamp').
    addEventListener("change", (e) => { showInSAMP(e.target.value); });

  document.getElementById('showinjs9').
    addEventListener("change", (e) => { showInJS9(e.target.value); });

  // Set up the button handlers, if we have any
  //
  if (settings.npages > 1) {
    for (let i = 1; i <= settings.npages; i++) {
      document.getElementById("page" + i.toString())
        .addEventListener("change", (e) => { setPage(e.target.value); });
    }
  }

  setPage(1);

  // Show/hide the color mapping table.
  //
  const colorButton = document.getElementById('colormappingselect');
  if (colorButton !== null) {
    colorButton
      .addEventListener("click",
			(e) => { toggleColorMapping(); });
  }

  // Update the "user content" section
  const usernotes = document.getElementById("usercontent");

  // We always use the user version if set, even if it is empty.
  //
  if (state.usernotes.user !== null) {
    usernotes.value = state.usernotes.user;
  } else {
    usernotes.value = state.usernotes.proposed;
  }

  // Do we allow the user to change this?
  // TODO: should also worry about the revision value
  //
  const saveusercontent = document.getElementById("saveusercontent");
  if (state.ensemble_status === "todo") {
      saveusercontent
	  .addEventListener("click",
			    (e) => { saveUserContent(); });
  } else {
      saveusercontent.disabled = true;
      document.getElementById("usercontent")
          .disabled = true;
  }

  // Now for the user action button.
  //
  // The choice is only present for the latest version, otherwise
  // it's just a label (i.e. non-interactive part of the UI).
  //
  let ua = "";
  if (state.useraction.user === null) {
    ua = state.useraction.proposed;
  } else {
    ua = state.useraction.user;
  }

  const choice = document.getElementById('choice');
  if (choice === null) {
      let lbl;
      if (ua === "accept") { lbl = "Accept"; }
      else if (ua === "delete") { lbl = "Delete Master"; }
      else if (ua === "manual") { lbl = "Manual"; }
      else if (ua === "") { lbl = "ERROR: no action recorded"; }
      else { lbl = ua; }

      document.getElementById('chosen').innerHTML = lbl;
      return;
  }

  // set up the chosen value
  choice.value = ua;

  if (state.ensemble_status === "todo") {
      choice.addEventListener("change",
			      (e) => { saveUserChoice(e.target.value); });
  } else {
      choice.disabled = true;
  }

}


const spinopts = {
    color: "#FF0000", opacity: 0.4,
    radius: 20, length: 20, width: 7
};


function initialize(opts) {

  settings = Object.assign({}, opts);

  /*** For the moment comment out the samp code, as the
       continual checking is making some unrelated
       tasks hard to debug/work on.
  if (typeof sampConnection === "undefined") {
    sampConnection = new samp.Connector("CHS Sender");
    sampConnection.onHubAvailability(sampIsAvailable, 2000);
  }
  ***/
  sampIsAvailable(false);

  let httpRequest = new XMLHttpRequest();
  if (!httpRequest) {
      alert("Unable to create a XMLHttpRequest!");
      return;
  }

  // Add the spinner to the whole page
  //
  let body = document.getElementsByTagName("body")[0];
  let spinner = new Spinner(spinopts);
  spinner.spin(body);

  httpRequest.addEventListener("load", function() {
    updatePage(httpRequest.response);
  });
  httpRequest.addEventListener("error", function() {
      alert("Unable to load data!");
  });
  httpRequest.addEventListener("abort", function() {
      alert("Unable to load data!");
  });
  httpRequest.addEventListener("loadend", function() {
      spinner.stop();
  });

  // Do I need to add a cache-busting identifier?
  httpRequest.open('GET', '/api/' +
		   settings.ensemble + '/' +
		   settings.revstr + '/' +
		   settings.masterid + '?' +
                   (new Date()).getTime());
  httpRequest.responseType = 'json';
  httpRequest.send();

}
