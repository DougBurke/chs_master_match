/*
 * Code for the master-hull page.
 */

"use strict";

var settings;

var state;

var pageNum;
var imageScale = 'log10';

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

// Store the region information for this *ensemble*; that is,
// provide coordinates (WCS) for all master hulls and stacks.
//
// This allows the UI to then decide whether to show all this
// information or not.
//
/***
var regionstore = {
  'masterhulls': {{ hull_store|safe }},
  'stackhulls': {{ stack_polys|safe }}
};
***/

// pan to this position; opts is the argument to pass
// to JS9 commands to determine the window to use.
// wcs contains ra and dec fields.
//
function goToRaDec(wcs, opts) {
  const pix = JS9.WCSToPix(wcs.ra, wcs.dec, opts);
  JS9.SetPan(pix.x, pix.y, opts);
}

// Finalize the widgets in the window and add the regions.
//
function finalizeJS9Display(stack, cptnum, winid) {

  return function(img) {
    setupJS9(img, stack, cptnum, winid);
    addRegionsToJS9(img, stack, winid);

    /* TODO: updateJS9StackCounter(winid, out.nhulls); */
/***
    const out = addRegionToJS9(img, ensemble, masterid, stack, winid,
                             data);
    updateJS9StackCounter(winid, out.nhulls);
    const opts = {display: img};
    if (typeof(out.center) !== "undefined") {
      goToRaDec(out.center, opts);
    }
***/
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
  //   - reload region
  //   - zoom
  //
  const opts = {display: img};
  let btns = container.getElementsByClassName("sigma");
  for (let i = 0; i < btns.length; i++) {
    btns[i].addEventListener("change", (e) => {
      JS9.GaussBlurData(e.target.value, opts);
    });
  }

  const sigma0 = document.getElementById(winid + "sigma0");

  btns = container.getElementsByClassName("binsize");
  for (let i = 0; i < btns.length; i++) {
    btns[i].addEventListener("change", (e) => {

      // What is the current position - I am not too bothered
      // if I do not calculate the exact center, as there is
      // an effective zoom in/out going on here due to the
      // change in binning - but let's see if this is visually
      // confusing (the issue is if I am really calculating the
      // center of the image here).
      //
      const idata = JS9.GetImageData(false, opts);

      // This does not seem to be working as intended!
      const wcs = JS9.PixToWCS(idata.width / 2, idata.height / 2, opts);

      JS9.DisplaySection({bin: e.target.value}, opts);

      // Could perhaps just jump to the master hull
      // goToRaDec(wcs, opts);  not working well

      // reset the blur button to 0 since rebinning removes
      // the blurring automatically
      sigma0.checked = true;
    });
  }

  // When reloading the regions, delete all the existing ones
  // as a precaution. The region info is reloaded from disk
  // (an alternative would be to just cache the data and re-use
  // it).
  //
  const reloadButton = document.getElementById(winid + 'ReloadRegions');
  if (reloadButton !== null) {
      // TODO: fix this (i.e. to match the latest behavior)
      reloadButton
	  .addEventListener("click", (e) => {
		  JS9.RemoveRegions("all", opts);
		  $.ajax({url: '/regions/ensemble/' +
			      settings.ensemble + '/' +
			      settings.masterid.toString(),
			      dataType: 'json'})
		  .done((data, textStatus) => {
			  addRegionToJS9(img, stack, winid, data);
		      });
	      });
  }

  // This currently doesn't work very well, as it loses information
  // on any other bin/... operation that has been applied.
  //
  const bandSelect = document.getElementById(winid + 'BandChoice');
  if (bandSelect !== null) {
      bandSelect
	  .addEventListener("change", (e) => {
		  const enfilter = band_to_filter(e.target.value);
		  console.log("new filter = [" + enfilter + "]");
		  JS9.DisplaySection({filter: enfilter}, opts);
	      });
  }

  const toggleButton = document.getElementById(winid + 'TogglePSFs');
  if (toggleButton !== null) {
      toggleButton
	  .addEventListener("click", (e) => { togglePSFs(stack, winid); });
  }

  const psfColSelect = document.getElementById(winid + 'PSFColor');
  if (psfColSelect !== null) {
      psfColSelect
	  .addEventListener("change", (e) => {
		  colorizePSFs(winid, e.target.value);
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
      .addEventListener("click", (e) => { JS9.DisplayPlugin('panner',
                                                            {display: winid}); });

  // Set up the layers. We have to base the options on
  // JS9.Regions.opts to get "sensible" behavior.
  //
  let layerOpts = Object.assign({}, JS9.Regions.opts);
  for (let name of ['movable', 'rotatable', 'resizable', 'evented']) {
    layerOpts[name] = false;
  }

  // Make sure we remove the onchange handler from these extra layers
  layerOpts.onchange = null;

  // Eric hopes this will address the half-pixel offset on rebin
  // but I need to investigate it further to find out what is
  // going in.
  //
  // layerOpts.dowcsstr = true;

  for (let name of [psfLayer, stackLayer, originalLayer, convexLayer]) {
    JS9.NewShapeLayer(name, layerOpts, opts);
  }
}

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
  JS9.AddShapes(layer, shape, opts, win);
  // JS9.AddRegions(shape, opts, win);

}


// Are the PSFs being shown or hidden for a stack?
var psfState = {};

function addPSFRegions(stack, win) {
  const psfs = settings.regionstore.stackpsfs[stack];
  if (typeof psfs === "undefined") {
    return;
  }

  const psfOpts = {color: 'yellow',
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

  psfState[stack] = true;
}

function togglePSFs(stack, winid) {
  const flag = !psfState[stack];
  JS9.ShowShapeLayer(psfLayer, flag, {display: winid});
  psfState[stack] = flag;
  let label = " PSFs";
  if (flag) { label = "Hide" + label; } else { label = "Show" + label; }
  document.getElementById(winid + 'TogglePSFs').innerHTML = label;
}

function colorizePSFs(winid, newcol) {
  JS9.ChangeShapes(psfLayer, "all", {color: newcol}, {display: winid});
}

// Add the stack-level and master hull(s) to the JS9 window.
//
// The PSF regions are drawn first (if available).
//
// Stack level hulls are drawn first (maybe in a different layer?)
//
function addRegionsToJS9(img, stack, regions) {

  let stackhulls = settings.regionstore.stackhulls[stack];
  if (typeof stackhulls === "undefined") {
    alert("no stack-level hulls found for " + stack);
    return;
  }

  let display = {display: img};

  addPSFRegions(stack, display);

  let linestyle = [1];
  let hullOpts = {color: 'orange',
                  strokeDashArray: linestyle,
                  changeable: false,
                  tags: 'stack'};

  for (let shull of stackhulls) {
    // linestyle: solid for mancode is 0, otherwise
    // dotted.
    //
    // wanted to add strokeWidth to increase the width, but can
    // not get it to work
    //
    if (shull.mancode === 0) {
      delete hullOpts.strokeDashArray;
    } else {
      hullOpts.strokeDashArray = linestyle;
    }
    add_hull_to_js9(shull, hullOpts, display, stackLayer);
  }

  /*
   * TODO: the tag should act as a unique identifier
   *       which will be useful when supporting QA cases
   *       (when can have multiple hulls).
   */
  let masterhull = settings.regionstore.masterhulls[settings.masterid];
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

  hullOpts = {movable: false,
              rotatable: false,
              resizable: false,
              tags: 'master'};

  /* If this is a review page then don't let the user change the hull */
  if (state.ensemble_status !== "todo") {
    hullOpts.changeable = false;
  }

  let ras = [];
  let decs = []

  const origOpts = Object.assign({}, hullOpts);
  origOpts.color = 'white';
  origOpts.strokeDashArray = [3, 3];

  for (let hull of masterhull.wcs) {
    add_hull_to_js9(hull, origOpts, display, originalLayer);
    add_hull_to_js9(hull, hullOpts, display, masterLayer);
    ras.push(hull.ra0);
    decs.push(hull.dec0);
  }

  // Pick the "middle" point if there are multiple master hulls (e.g. QA).
  // This fails if the hulls straddle ra=0/360.
  //
  let ra0 = 0;
  let dec0 = 0;
  if (ras.length > 1) {
    ra0 = 0.5 + (Math.min(...ras) + Math.max(...ras));
    dec0 = 0.5 + (Math.min(...decs) + Math.max(...decs));
  } else {
    ra0 = ras[0];
    dec0 = decs[0];
  }
  goToRaDec({ra: ra0, dec: dec0}, display);
}

// The following is OLD
//
// Returns information useful when the display is being created,
// but not when the regions are being re-created.
//
// TODO: rework this now moving knowledge into JS
//
function addRegionToJS9(img, stack, winid, regions) {

  /* can have multiple components for a stack */
  let linestyle;
  let nstackhulls = 0;

  const js9win = {display: img};

  const hulls = [];
  for (let i = 0; i < regions.stacks.length; i++) {
    let shull = regions.stacks[i];
    if (shull.stack === stack) {

      hulls.push(shull);

      /* experiment with making stack-level hulls "fixed", but
         I can imagine this could be annoying at times (useful
         at others) */
      linestyle = [1];  /* how to change this */

      add_hull_to_js9(shull,
                      {color: 'orange',
                       strokeDashArray: linestyle,
                       changeable: false,
                       tags: 'stack'},
                      js9win);

      nstackhulls += 1;
    }
  }

  // Ensure that the store is updated
  if (hulls.length > 0) {
    settings.regionstore[stack] = hulls;
  } else {
    delete settings.regionstore[stack];
  }

  /*
   * Add master hull last (so it appears on top).
   * Trying to add the right level of edit-ability to the polygon.
   *
   * Note that points can be moved and deleted even with these options
   * set.
   */

  if (typeof regions.master !== "undefined") {

    add_hull_to_js9(regions.master,
                    {movable: false,
                     rotatable: false,
                     resizable: false,  // do we want this true?
                     tags: 'master'},
                    js9win);

    settings.regionstore.master = [regions.master];
  } else {
    delete settings.regionstore.master;
  }

  const out = {nhulls: nstackhulls};
  if (typeof(regions.center) !== "undefined") {
    out.center = regions.center;
  }
  return out;
}

// Add the number of stacks to the JS9 window.
function updateJS9StackCounter(winid, n) {
  document.getElementById(winid + 'stackid')
    .innerHTML += "<span style='float: right;'>Stack hulls: " +
                  n.toString() + "</span>";
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
      "</span></div>" +
    "<div class='useropts'>";

  html += "<div class='buttonopts'>";
  let name = id + "sigma";
  html += "<div class='blur'><span>Blur</span>";
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

  const all_bins = [1, 2, 4, 8, 16, 32, 64, 128];
  let def_binsize, def_start;
  if (stack.startsWith('hrc')) {
    def_binsize = 64;
    def_start = 1;
  } else {
    def_binsize = 8;
    def_start = 0;
  }
  const def_bins = all_bins.slice(def_start, def_start + 7);

  name = id + "binsize";
  html += "<div class='rebin'><span>Bin</span>";
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
  html += "</div>";  // class=buttonopts

  /* hide for now
  html += "<div class='reload'>";
  html += "<button id='" + id;
  html += "ReloadRegions'>Reload Regions</button>";
  html += "</div>";
  */

  /* hide until the behavior of DisplaySection has either been improved
     or a way to get at the other information needed has been added to
     the JS9 API
  if (band !== "w") {
      html += "<div class='bandchoice'>Band: ";
      html += "<select id='" + id + "BandChoice'>";

      for (const bandval of ['b', 'u', 's', 'm', 'h']) {
	  html += "<option value='" + bandval + "'";
	  if (bandval === band) {
	      html += " selected";
	  }
	  html += ">" + bandval + "</option>";
      }

      html += "</select>";
      html += "</div>";
  }
  */

  const psfs = settings.regionstore.stackpsfs[stack];
  if (typeof psfs !== "undefined") {
      html += "<div class='reload'>";
      html += "<button id='" + id;
      html += "TogglePSFs'>Hide PSFs</button>";
      html += "</div>";

      html += "<div class='colorize'>";
      html += "<select id='" + id + "PSFColor'>";
      for (const newcol of ['yellow', 'white', 'black', 'red',
			    'orange', 'cyan', 'blue', 'brown']) {
	  html += "<option value='" + newcol + "'";
	  if (newcol === 'yellow') {
	      html += " selected";
	  }
	  html += ">" + newcol + "</option>";
      }

      html += "</select>";
      html += "</div>";

  }


  html += "<div class='zoom'>Zoom: ";
  html += "<button class='zoomin'>In</button>";
  html += "<button class='zoomout'>Out</button>";
  html += "</div>";

  html += "<div class='showable'>"; // "Toggle: ";
  html += "<button id='" + id + "ShowPanner'>Panner</button>";
  html += "</div>";

  html += "</div>";

  html += "<div class='JS9Menubar' id='" + id + "Menubar'></div>";
  html += "<div class='JS9' id='" + id + "'></div>'";
  html += "<div class='JS9Colorbar' id='" + id + "'></div>'";

  return html;
}


// "unique" id for new JS9 windows
//
//
var idctr = 1;
function js9_id(stack, cptnum) {
  // return stack;

  var retval = idctr.toString();
  idctr += 1;
  return "js9win" + retval;
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
    opts.bin = 8;
    opts.xdim = 8192;
    opts.ydim = 8192;
    // does setting this actually do anything?
    opts.filter = band_to_filter(band);
  } else {
    opts.bin = 64;
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
  if (action.tags[0] !== "master") { return; }

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

  // Convert the polygon into a convex hull (if necessary).
  //
  let chull_sky = window.convexHull(action.pts);
  let chull_eqpos = [];
  for (let sky of chull_sky) {
    chull_eqpos.push(JS9.PixToWCS(sky.x, sky.y, baseWin));
  }

  // Now that the convex-hull version appears in a different layer,
  // the name can match the master name.
  //
  // const convexName = 'convex';
  const convexName = 'master';

  const convexOpts = {color: 'cyan',
                      strokeDashArray: [5, 3],
                      changeable: false,
                      evented: false,
                      tags: convexName};

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
    let hdl = JS9.GetImage(imname);

    // stop propogating the onchange signal
    hdl.params.xeqonchange = !hdl.params.xeqonchange;

    // Do we have a convex hull already?
    // let rs = JS9.GetRegions(convexName, imname);
    let rs = JS9.GetShapes(convexLayer, convexName, imname);
    if (rs.length === 0) {
      var ras = [];
      var decs = [];
      for (let wcs of chull_eqpos) {
        ras.push(wcs.ra);
        decs.push(wcs.dec);
      }
      add_hull_to_js9({ra: ras, dec: decs}, convexOpts,
                      imname, convexLayer);
    } else {
      const hullpts = [];
      for (let wcs of chull_eqpos) {
        hullpts.push(JS9.WCSToPix(wcs.ra, wcs.dec, imname));
      }

      // JS9.ChangeRegions(convexName, {pts: hullpts}, imname);
      JS9.ChangeShapes(convexLayer, convexName, {pts: hullpts},
                       imname);
    }

    // Only adjust the master polygon if this is not the
    // window the user is changing.
    if (owin !== js9win) {

      // Need to convert to image coordinates
      const pts = [];
      for (let wcs of action.wcspts) {
        pts.push(JS9.WCSToPix(wcs.ra, wcs.dec, imname));
      }

      // JS9.ChangeRegions('master', {pts: pts}, imname);
      JS9.ChangeShapes(masterLayer, 'master', {pts: pts}, imname);
    }

    hdl.params.xeqonchange = !hdl.params.xeqonchange;

    // Ensure the regions layer is where the user action happens
    JS9.ActiveShapeLayer(masterLayer, imname);
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

// Set up the page

function updatePage(json) {

  state = Object.assign({}, json);

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
