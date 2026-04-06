// Shawn Gano - Stellarium Functions 
// 17 Feb 2026
//   place in same folder as Orion_ArtemisII.ssc script
//


//------------------------------------------------------------------
// Convert Declination   degrees to string (e.g., "+45d45m30.23s")
function degToDmsString(degDecimal) 
{
    // store decimal degrees
	let deg =  degDecimal
    
    // Handle negative angles
    let sign = deg < 0 ? "-" : "+";
    deg = Math.abs(deg);
    
    // Calculate Degrees, Minutes, Seconds
    let d = Math.floor(deg);
    let m = Math.floor((deg - d) * 60);
    let s = ((deg - d - m / 60) * 3600); // 2 decimal places for seconds  // moved .toFixed(2) to string versions

    // String versions with padding
    let dStr  = String(  d  ).padStart(2, '0'); // not ensures at least two digist, if number is more than two no change
    let mStr = String ( m ).padStart(2, '0');
    let sStr  =  String(  s.toFixed(2)  ).padStart(5, '0'); // 5 digits since: 00.00 [decimal is one char]

    // Format string
    //return sign + d + "d" + m + "m" + s + "s";
    return sign + dStr + "d" + mStr + "m" + sStr + "s";

} // degToDmsString --------------------------
// Example usage in Stellarium script:
//var degreesDec = 1.5708; 
//core.debug( degToDmsString(degreesDec) );  // output:  +1d34m14.88s
//var degreesDec2 = 101.12345; 
//core.debug( degToDmsString(degreesDec2) );  // output:  +101d07m24.42s (invalid value- but testing above 100)
//var degreesDec2 = -20.0022; 
//core.debug( degToDmsString(degreesDec2) );  // output:  -20d00m07.92s
//------------------------------------------------------------------


//------------------------------------------------------------------
// Convert Right Ascention hours (double) to HMS String (e.g., "01h06m39.96s" )
function raHoursToHmsString(raHours) 
{
    // Ensure hours are within 0-24 range (well at least non-negative)
    let h = raHours % 24;
    if (h < 0) h += 24;
    
    // Extract hours, minutes, seconds
    let hours = Math.floor( h );
    let minutes = Math.floor( (h - hours) * 60 );
    let seconds = ((h - hours) * 60 - minutes) * 60;

    // string versionsof fixed length
    let hoursStr     = String(  hours  ).padStart(2, '0');
    let minutesStr = String(  minutes  ).padStart(2, '0');
    let secondsStr = String(   seconds.toFixed(2)  ).padStart(5, '0'); // 5 digits since: 00.00 [decimal is one char]
    
    // Format as "HHhMMmSSs" (padded to 2 digits) ("11h19m07.7s")
    //return `${String(hours).padStart(2, '0')}h${String(minutes).padStart(2, '0')}m${String(Math.round(seconds)).padStart(2, '0')}s`;
    return   hoursStr + "h" + minutesStr + "m" + secondsStr + "s";

} // raHoursToHMS
// Example usage:
//var ra = 10.5; // hours
//core.debug(raHoursToHMS(ra)); // Outputs: 10h30m00.00s
//var ra2 = 1.1111; // hours
//core.debug(raHoursToHMS(ra2)); // Outputs: 01h06m39.96s
//------------------------------------------------------------------


/**
 * Performs linear interpolation within a sorted lookup table (array of [x, y] points).
 *
 * @param {Array<Array<number>>} data - Sorted array of data points, e.g., [[x0, y0], [x1, y1], ...].
 * @param {number} x - The x value to interpolate the corresponding y value for.
 * @returns {number | null} The interpolated y value, or null if x is out of the table's range.
 */
function lookupTableLinearInterpolation(data, x) 
{
	// Check if x is within the bounds of the table
  	if (x < data[0][0] || x > data[data.length - 1][0]) 
	{
    		// For values outside the defined range, one could also extrapolate or clamp as needed
    		return null;
  	}
	
  	// Find the two data points that surround the given x value
  	for (let i = 0; i < data.length - 1; i++) 
	{
    		const x0 = data[i][0];
    		const y0 = data[i][1];
    		const x1 = data[i + 1][0];
    		const y1 = data[i + 1][1];

		// check to see if these two surround given value
    		if (x0 <= x && x <= x1) 
		{
      			// Perform linear interpolation using the formula: y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
      			const t = (x - x0) / (x1 - x0); // The interpolation factor (0.0 to 1.0)
      			const interpolatedY = y0 + (y1 - y0) * t;

      			return interpolatedY;
    		} // if found the two bounding data points

  	} // for
  
  // Should be unreachable if the initial range check passes and data is sorted
  return null; 
} // lookupTableLinearInterpolation

// ---------------------------------------------------------
