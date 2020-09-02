from distutils.core import setup, Extension


setup(name='sample',
      ext_modules=[
        Extension('sharedmatting',
                  ['wrapper.cpp', 'sharedmatting.cpp'],
                  )
        ]
)
