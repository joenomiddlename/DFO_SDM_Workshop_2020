$(document).ready(function() {
  $plots = $('img.plot');
  $chunks = $('pre').has('code');
  $chunks = $chunks.filter(function(idx) {
    return $(this).children('code').outerHeight(false) > parseInt($(this).css('line-height'));
  });
  
  $chunks.each(function () {
    if($(this).hasClass('r')) {
      $(this).append("<div class=\"showopt\">Show Answer Code</div><br style=\"line-height:22px;\"/>");
    } else {
      $(this).append("<div class=\"showopt\">Show Answer Output</div><br style=\"line-height:22px;\"/>");
    }
  });
  
  $plots.each(function () {
    $(this).wrap('<pre class=\"plot\"></pre>');
    $(this).parent('pre.plot').prepend("<div class=\"showopt\">Show Answer Plot</div><br style=\"line-height:22px;\"/>");
  });
  
  // hide all chunks when document is loaded
  $chunks.children('code').toggle();
  $('pre.plot').children('img').toggle();
  // function to toggle the visibility
  $('.showopt').click(function() {
    var label = $(this).html();
    if (label.indexOf("Show") >= 0) {
      $(this).html(label.replace("Show", "Hide"));
    } else {
      $(this).html(label.replace("Hide", "Show"));
    }
    $(this).siblings('code, img').slideToggle('fast', 'swing');
  });
});