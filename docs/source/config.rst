Pycoast from a configuration file
---------------------------------

If you want to run to avoid typing the same options over and over again, or if
caching is an optimization you want, you can use a configuration file with the
pycoast options you need:

.. code-block:: ini

  [cache]
  file=/var/run/satellit/white_overlay
  regenerate=False

  [coasts]
  level=1
  width=0.75
  outline=white
  fill=yellow

  [borders]
  outline=white
  width=0.5

Then, you can just call:

   >>> cw.add_overlay_from_config("my_config.cfg", area_def)
